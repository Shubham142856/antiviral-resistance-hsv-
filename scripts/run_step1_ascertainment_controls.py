"""
run_step1_ascertainment_controls.py
Executes Step 1 Ascertainment Bias & Curation Controls under the pre-registered protocol:
docs/STEP1_ASCERTAINMENT_CONTROL_PREREGISTRATION.md (Commit 756ffe6, SHA-256 298c59cbe...)
"""
import re
import json
import pathlib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# 1. Load data
feat_df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
herpes = pd.read_csv('herpesdrg-db/herpesdrg-db.tsv', sep='\t')
hsv1 = herpes[herpes['virus'] == 'HSV1'].copy()

# Parse aa_change in HerpesDRG
def parse_mut(s):
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', str(s).strip())
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    return None, None, None

parsed = [parse_mut(x) for x in hsv1['aa_change']]
hsv1['wt'] = [p[0] for p in parsed]
hsv1['pos'] = [p[1] for p in parsed]
hsv1['mut'] = [p[2] for p in parsed]
hsv1['gene_std'] = hsv1['gene'].replace({'TK': 'UL23', 'POL': 'UL30', 'DPOL': 'UL30'})

# Deduplicate herpesdrg per (gene_std, pos, wt, mut) prioritizing recombinant/marker transfer if present
def get_tier_rank(sg):
    s = str(sg).lower()
    if 'marker transfer' in s or 'recombinant' in s or 'lab derived' in s or 'lab induced' in s:
        return 2
    if 'isolated' in s:
        return 1
    return 0

hsv1['tier_rank'] = [get_tier_rank(x) for x in hsv1['strain_generation']]
hsv1_sorted = hsv1.sort_values(by=['gene_std', 'pos', 'wt', 'mut', 'tier_rank'], ascending=[True, True, True, True, False])
hsv1_dedup = hsv1_sorted.drop_duplicates(subset=['gene_std', 'pos', 'wt', 'mut']).copy()

# Merge onto feat_df
merged_df = pd.merge(
    feat_df,
    hsv1_dedup[['gene_std', 'pos', 'wt', 'mut', 'strain_generation', 'test_method', 'Aciclovir', 'Foscarnet', 'ref_doi']],
    left_on=['gene', 'pos', 'wt', 'mut'],
    right_on=['gene_std', 'pos', 'wt', 'mut'],
    how='left'
)

# Spec-exact Where
for g in ['UL23', 'UL30']:
    m = (merged_df['gene'] == g)
    z_mu = - (merged_df.loc[m, 'mu_plm_19'] - merged_df.loc[m, 'mu_plm_19'].mean()) / merged_df.loc[m, 'mu_plm_19'].std()
    z_dist = (merged_df.loc[m, 'dist_to_site_angstrom'] - merged_df.loc[m, 'dist_to_site_angstrom'].mean()) / merged_df.loc[m, 'dist_to_site_angstrom'].std()
    merged_df.loc[m, 'where_specexact'] = (z_mu - z_dist) / np.std(z_mu - z_dist)
    merged_df.loc[m, 'constraint_only'] = z_mu
    merged_df.loc[m, 'proximity_only'] = - z_dist

def is_numeric_fc(val):
    if pd.isna(val): return False
    s = str(val).replace('>', '').replace('<', '').strip()
    try:
        float(s)
        return True
    except:
        return False

def get_evidence_tier(sg):
    if pd.isna(sg): return 'Unspecified / Literature'
    s = str(sg).strip().lower()
    if 'marker transfer' in s or 'recombinant' in s or 'lab induced' in s or 'lab derived' in s:
        return 'Recombinant / Lab-Validated'
    if 'isolated strain' in s:
        return 'Clinical Isolate'
    return 'Other'

merged_df['evidence_tier'] = [get_evidence_tier(x) for x in merged_df['strain_generation']]

rng = np.random.RandomState(42)

def compute_stratum_metrics(sub_df, score_col='where_specexact', n_boot=1000):
    y = sub_df['target_binary'].values.astype(int)
    scores = sub_df[score_col].values
    pos = sub_df['pos'].values
    
    n_total = len(y)
    n_pos = int(y.sum())
    n_neg = n_total - n_pos
    
    if n_pos == 0 or n_neg == 0:
        return {'n': n_total, 'pos': n_pos, 'neg': n_neg, 'auroc': None, 'ci_95': None, 'status': 'Undefined (single class)'}
        
    auc = float(roc_auc_score(y, scores))
    
    # Position-clustered bootstrap
    u_pos = np.unique(pos)
    pos_map = {p: np.where(pos == p)[0] for p in u_pos}
    boot_aucs = []
    for _ in range(n_boot):
        s_p = rng.choice(u_pos, size=len(u_pos), replace=True)
        idx = np.concatenate([pos_map[p] for p in s_p])
        y_b = y[idx]
        if len(np.unique(y_b)) > 1:
            boot_aucs.append(float(roc_auc_score(y_b, scores[idx])))
            
    ci = [float(np.percentile(boot_aucs, 2.5)), float(np.percentile(boot_aucs, 97.5))] if len(boot_aucs) > 50 else [None, None]
    
    return {
        'n': n_total,
        'pos': n_pos,
        'neg': n_neg,
        'unique_positions': len(u_pos),
        'auroc': auc,
        'ci_95': ci
    }

results = {
    'protocol': {
        'document': 'docs/STEP1_ASCERTAINMENT_CONTROL_PREREGISTRATION.md',
        'commit': '756ffe6',
        'sha256': '298c59cbe64e06aece43d179287a55190f70e34d9a36194ed9d58adadca22937'
    },
    'targets': {}
}

targets = [
    ('UL30', 'Aciclovir', 'Primary Target'),
    ('UL30', 'Foscarnet', 'Secondary Target'),
    ('UL23', 'Aciclovir', 'Secondary Target')
]

print("="*105)
print("STEP 1: ASCERTAINMENT BIAS & CURATION CONTROL EVALUATION")
print("="*105)

for g, d, role in targets:
    sub = merged_df[(merged_df['gene'] == g) & (merged_df['drug'] == d)].copy().reset_index(drop=True)
    full_metrics = compute_stratum_metrics(sub, 'where_specexact')
    full_auc = full_metrics['auroc']
    
    target_res = {
        'role': role,
        'full_cohort': {
            'where': full_metrics,
            'constraint_only': compute_stratum_metrics(sub, 'constraint_only'),
            'proximity_only': compute_stratum_metrics(sub, 'proximity_only')
        },
        'strata': {}
    }
    
    print(f"\nTarget: {g}-{d} ({role}) | Full Cohort AUROC = {full_auc:.4f} (N={full_metrics['n']}, Pos={full_metrics['pos']})")
    print("-" * 105)
    print(f"{'Stratum / Subset':<32} | {'N (Pos/Neg)':<14} | {'Where AUROC [95% CI]':<26} | {'Constraint':<10} | {'Proximity':<10} | {'Decision'}")
    print("-" * 105)
    
    # 1. Measured fold-change subset
    sub['has_num_fc'] = [is_numeric_fc(x) for x in sub[d]]
    num_sub = sub[sub['has_num_fc']].copy().reset_index(drop=True)
    m_num = compute_stratum_metrics(num_sub, 'where_specexact')
    m_num_c = compute_stratum_metrics(num_sub, 'constraint_only')
    m_num_p = compute_stratum_metrics(num_sub, 'proximity_only')
    
    if m_num['auroc'] is not None and m_num['pos'] >= 10:
        delta = m_num['auroc'] - full_auc
        decision = "PASS" if delta >= -0.05 else ("INCONCLUSIVE" if delta >= -0.10 else "ARTIFACT")
    else:
        decision = "UNDERPOWERED (<10 events)"
        
    target_res['strata']['measured_fold_change_subset'] = {
        'where': m_num,
        'constraint_only': m_num_c,
        'proximity_only': m_num_p,
        'delta_vs_full': (m_num['auroc'] - full_auc) if m_num['auroc'] is not None else None,
        'decision': decision
    }
    ci_str = f"[{m_num['ci_95'][0]:.4f}, {m_num['ci_95'][1]:.4f}]" if m_num['ci_95'][0] is not None else "N/A"
    n_str = f"{m_num['n']} ({m_num['pos']}/{m_num['neg']})"
    auc_str = f"{m_num['auroc']:.4f} {ci_str}"
    c_str = f"{m_num_c['auroc']:.4f}"
    p_str = f"{m_num_p['auroc']:.4f}"
    print(f"{'Measured Fold-Change Only':<32} | {n_str:<14} | {auc_str:<26} | {c_str:<10} | {p_str:<10} | {decision}")
    
    # 2. Evidence Tiers
    for tier in ['Recombinant / Lab-Validated', 'Clinical Isolate']:
        tier_sub = sub[sub['evidence_tier'] == tier].copy().reset_index(drop=True)
        m_t = compute_stratum_metrics(tier_sub, 'where_specexact')
        m_tc = compute_stratum_metrics(tier_sub, 'constraint_only')
        m_tp = compute_stratum_metrics(tier_sub, 'proximity_only')
        
        if m_t['auroc'] is not None and m_t['pos'] >= 10:
            delta = m_t['auroc'] - full_auc
            decision = "PASS" if delta >= -0.05 else ("INCONCLUSIVE" if delta >= -0.10 else "ARTIFACT")
        else:
            decision = f"UNDERPOWERED ({m_t['pos']} events)"
            
        target_res['strata'][tier] = {
            'where': m_t,
            'constraint_only': m_tc,
            'proximity_only': m_tp,
            'delta_vs_full': (m_t['auroc'] - full_auc) if m_t['auroc'] is not None else None,
            'decision': decision
        }
        ci_str = f"[{m_t['ci_95'][0]:.4f}, {m_t['ci_95'][1]:.4f}]" if m_t['ci_95'] and m_t['ci_95'][0] is not None else "N/A"
        n_str = f"{m_t['n']} ({m_t['pos']}/{m_t['neg']})"
        auc_str = f"{m_t['auroc']:.4f} {ci_str}" if m_t['auroc'] is not None else "N/A"
        c_str = f"{m_tc['auroc']:.4f}" if m_tc['auroc'] is not None else "N/A"
        p_str = f"{m_tp['auroc']:.4f}" if m_tp['auroc'] is not None else "N/A"
        print(f"{tier:<32} | {n_str:<14} | {auc_str:<26} | {c_str:<10} | {p_str:<10} | {decision}")
        
    # 3. Assay methods (PRA vs Dye Uptake)
    for method in ['PRA', 'dye uptake assay']:
        meth_sub = sub[sub['test_method'] == method].copy().reset_index(drop=True)
        m_m = compute_stratum_metrics(meth_sub, 'where_specexact')
        m_mc = compute_stratum_metrics(meth_sub, 'constraint_only')
        m_mp = compute_stratum_metrics(meth_sub, 'proximity_only')
        
        if m_m['auroc'] is not None and m_m['pos'] >= 10:
            delta = m_m['auroc'] - full_auc
            decision = "PASS" if delta >= -0.05 else ("INCONCLUSIVE" if delta >= -0.10 else "ARTIFACT")
        else:
            decision = f"UNDERPOWERED ({m_m['pos']} events)"
            
        target_res['strata'][method] = {
            'where': m_m,
            'constraint_only': m_mc,
            'proximity_only': m_mp,
            'delta_vs_full': (m_m['auroc'] - full_auc) if m_m['auroc'] is not None else None,
            'decision': decision
        }
        ci_str = f"[{m_m['ci_95'][0]:.4f}, {m_m['ci_95'][1]:.4f}]" if m_m['ci_95'] and m_m['ci_95'][0] is not None else "N/A"
        n_str = f"{m_m['n']} ({m_m['pos']}/{m_m['neg']})"
        auc_str = f"{m_m['auroc']:.4f} {ci_str}" if m_m['auroc'] is not None else "N/A"
        c_str = f"{m_mc['auroc']:.4f}" if m_mc['auroc'] is not None else "N/A"
        p_str = f"{m_mp['auroc']:.4f}" if m_mp['auroc'] is not None else "N/A"
        print(f"{method:<32} | {n_str:<14} | {auc_str:<26} | {c_str:<10} | {p_str:<10} | {decision}")
        
    results['targets'][f"{g}-{d}"] = target_res

pathlib.Path('logs').mkdir(exist_ok=True)
with open('logs/step1_ascertainment_controls.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "="*105)
print("Saved Step 1 results to logs/step1_ascertainment_controls.json")
