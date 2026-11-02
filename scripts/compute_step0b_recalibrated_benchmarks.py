"""
compute_step0b_recalibrated_benchmarks.py
Implements Step 0b refinements:
1. Strict pre-registered N_pos >= 5 rule for domain evaluability (FOS drops Exo [N_pos=2] and Thumb [N_pos=3]).
2. Excess AUROC reference: max(0.50, AUROC_pos_only).
3. Exact distinct domain-bootstrap resamples count (small-n proof).
4. Sign-consistency across evaluable domains (k/n positive).
5. Exact BLOSUM62 domain scores on UL23 (verifying domain-level truth vs legacy 5-block copy).
6. Bimodality diagnostic on beta5 (leave-domain-out influence).
7. Variant-drug records vs unique substitutions and position overlap.
Saves results to logs/step0b_recalibrated_benchmarks.json.
"""
import json
import math
import pathlib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

# 1. Load data
df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
spatial_df = pd.read_csv('data/processed/primary_3d_spatial_folds.tsv', sep='\t')
cols_to_merge = ['coord_x', 'coord_y', 'coord_z', 'has_3d_coord', 'score_esm2_650m', 'score_protbert', 'score_blosum62']
df = pd.merge(df, spatial_df[['gene', 'pos', 'wt', 'mut', 'drug'] + cols_to_merge], on=['gene', 'pos', 'wt', 'mut', 'drug'], how='left')

# Spec-exact Where
for g in ['UL23', 'UL30']:
    m = (df['gene'] == g)
    z_mu = - (df.loc[m, 'mu_plm_19'] - df.loc[m, 'mu_plm_19'].mean()) / df.loc[m, 'mu_plm_19'].std()
    z_dist = (df.loc[m, 'dist_to_site_angstrom'] - df.loc[m, 'dist_to_site_angstrom'].mean()) / df.loc[m, 'dist_to_site_angstrom'].std()
    df.loc[m, 'where_specexact'] = (z_mu - z_dist) / np.std(z_mu - z_dist)

def get_dual_buffer_mask(sub_df, test_domain, seq_buf=5, spatial_buf=12.0):
    te_mask = (sub_df['domain_5'].values == test_domain)
    test_positions = sub_df.loc[te_mask, 'pos'].values
    has_coord_test = sub_df[te_mask & (sub_df['has_3d_coord'] == 1)]
    tr_drop = np.zeros(len(sub_df), dtype=bool)
    for p in test_positions:
        tr_drop |= (np.abs(sub_df['pos'].values - p) <= seq_buf)
    if len(has_coord_test) > 0:
        tx = has_coord_test['coord_x'].values
        ty = has_coord_test['coord_y'].values
        tz = has_coord_test['coord_z'].values
        c_mask = (sub_df['has_3d_coord'] == 1).values
        cx = sub_df['coord_x'].values
        cy = sub_df['coord_y'].values
        cz = sub_df['coord_z'].values
        for i in range(len(has_coord_test)):
            dist = np.sqrt((cx - tx[i])**2 + (cy - ty[i])**2 + (cz - tz[i])**2)
            tr_drop |= (c_mask & (dist < spatial_buf))
    tr_mask = (~te_mask) & (~tr_drop)
    return tr_mask, te_mask

targets = [
    ('UL30', 'Aciclovir', 'Primary Target'),
    ('UL30', 'Foscarnet', 'Secondary Target (2 evaluable domains per N_pos>=5 rule)'),
    ('UL23', 'Aciclovir', 'Comparative Kinase (Dispensable Target)')
]

results = {}
rng = np.random.RandomState(42)

for g, d, desc in targets:
    sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
    domains = sorted(sub['domain_5'].unique())
    y_all = sub['target_binary'].values
    
    domain_records = []
    
    for dom in domains:
        tr_m, te_m = get_dual_buffer_mask(sub, dom)
        y_te = y_all[te_m]
        n_pos_test = int(y_te.sum())
        
        # STRICT PRE-REGISTERED RULE: N_pos >= 5 for evaluability
        if n_pos_test >= 5 and len(np.unique(y_te)) > 1 and len(np.unique(y_all[tr_m])) > 1:
            rf_pos = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=1)
            rf_pos.fit(sub.loc[tr_m, ['pos']], y_all[tr_m])
            p_pos = rf_pos.predict_proba(sub.loc[te_m, ['pos']])[:, 1]
            
            auc_p = float(roc_auc_score(y_te, p_pos))
            auc_esm = float(roc_auc_score(y_te, sub.loc[te_m, 'score_esm2_650m']))
            auc_bert = float(roc_auc_score(y_te, sub.loc[te_m, 'score_protbert']))
            auc_blo = float(roc_auc_score(y_te, sub.loc[te_m, 'score_blosum62']))
            auc_wh = float(roc_auc_score(y_te, sub.loc[te_m, 'where_specexact']))
            
            # REFINED REFERENCE: max(0.50, auc_p)
            ref_baseline = max(0.50, auc_p)
            
            domain_records.append({
                'domain': dom,
                'n_test': int(len(y_te)),
                'n_pos': n_pos_test,
                'auc_pos_rf': auc_p,
                'reference_baseline': ref_baseline,
                'auc_esm2': auc_esm,
                'auc_protbert': auc_bert,
                'auc_blosum62': auc_blo,
                'auc_where': auc_wh,
                'excess_esm2': auc_esm - ref_baseline,
                'excess_protbert': auc_bert - ref_baseline,
                'excess_blosum62': auc_blo - ref_baseline,
                'excess_where': auc_wh - ref_baseline
            })
            
    n_dom = len(domain_records)
    
    # Domain-clustered bootstrap with tracking of DISTINCT multisets
    boot_excess_esm = []
    boot_excess_wh = []
    boot_excess_bert = []
    boot_excess_blo = []
    distinct_multisets = set()
    
    for _ in range(2000):
        # Draw indices with replacement
        idx_resamp = tuple(sorted(rng.choice(range(n_dom), size=n_dom, replace=True)))
        distinct_multisets.add(idx_resamp)
        
        boot_excess_esm.append(np.mean([domain_records[i]['excess_esm2'] for i in idx_resamp]))
        boot_excess_wh.append(np.mean([domain_records[i]['excess_where'] for i in idx_resamp]))
        boot_excess_bert.append(np.mean([domain_records[i]['excess_protbert'] for i in idx_resamp]))
        boot_excess_blo.append(np.mean([domain_records[i]['excess_blosum62'] for i in idx_resamp]))
        
    ci_esm = [float(np.percentile(boot_excess_esm, 2.5)), float(np.percentile(boot_excess_esm, 97.5))]
    ci_wh = [float(np.percentile(boot_excess_wh, 2.5)), float(np.percentile(boot_excess_wh, 97.5))]
    ci_bert = [float(np.percentile(boot_excess_bert, 2.5)), float(np.percentile(boot_excess_bert, 97.5))]
    ci_blo = [float(np.percentile(boot_excess_blo, 2.5)), float(np.percentile(boot_excess_blo, 97.5))]
    
    mean_esm = float(np.mean([r['auc_esm2'] for r in domain_records]))
    mean_blo = float(np.mean([r['auc_blosum62'] for r in domain_records]))
    mean_wh = float(np.mean([r['auc_where'] for r in domain_records]))
    mean_pos = float(np.mean([r['auc_pos_rf'] for r in domain_records]))
    
    mean_excess_esm = float(np.mean([r['excess_esm2'] for r in domain_records]))
    mean_excess_blo = float(np.mean([r['excess_blosum62'] for r in domain_records]))
    mean_excess_wh = float(np.mean([r['excess_where'] for r in domain_records]))
    
    # Sign consistency
    sign_esm = f"{sum(r['excess_esm2'] > 0 for r in domain_records)}/{n_dom}"
    sign_blo = f"{sum(r['excess_blosum62'] > 0 for r in domain_records)}/{n_dom}"
    sign_wh = f"{sum(r['excess_where'] > 0 for r in domain_records)}/{n_dom}"
    
    results[f"{g}-{d}"] = {
        'description': desc,
        'n_evaluable_domains': n_dom,
        'evaluable_domains_list': [r['domain'] for r in domain_records],
        'distinct_bootstrap_multisets': len(distinct_multisets),
        'theoretical_multisets_max': int(math.comb(n_dom + n_dom - 1, n_dom)),
        'domain_records': domain_records,
        'domain_means': {
            'auc_pos_rf': mean_pos,
            'auc_esm2': mean_esm,
            'auc_blosum62': mean_blo,
            'auc_where': mean_wh
        },
        'excess_metrics': {
            'reference_rule': 'max(0.50, auc_pos_rf)',
            'esm2': {
                'domain_mean_excess': mean_excess_esm,
                'ci_95': ci_esm,
                'sign_consistency': sign_esm
            },
            'blosum62': {
                'domain_mean_excess': mean_excess_blo,
                'ci_95': ci_blo,
                'sign_consistency': sign_blo
            },
            'where': {
                'domain_mean_excess': mean_excess_wh,
                'ci_95': ci_wh,
                'sign_consistency': sign_wh
            }
        }
    }

# 2. Leave-One-Domain-Out Diagnostic on Beta5 (Bimodality / Influence Investigation)
ul30_sub = df[(df['gene'] == 'UL30') & (df['drug'] == 'Aciclovir')].copy().reset_index(drop=True)
W = ul30_sub['where_specexact'].values
S = - (ul30_sub['severity_plm'] - ul30_sub['severity_plm'].mean()) / ul30_sub['severity_plm'].std()
y = ul30_sub['target_binary'].values.astype(int)

W_std = (W - np.mean(W)) / np.std(W)
S_std = (S - np.mean(S)) / np.std(S)
X = np.column_stack([W_std, S_std, S_std**2, W_std * S_std, W_std * (S_std**2)])

loo_domain_beta5 = {}
for dom in sorted(ul30_sub['domain_5'].unique()):
    mask = (ul30_sub['domain_5'] != dom).values
    clf_loo = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42)
    clf_loo.fit(X[mask], y[mask])
    b5 = float(clf_loo.coef_[0][4])
    b1 = float(clf_loo.coef_[0][0])
    loo_domain_beta5[dom] = {'beta5': b5, 'beta1_where': b1, 'n_dropped': int((~mask).sum()), 'pos_dropped': int(y[~mask].sum())}

results['beta5_bimodality_diagnostic'] = {
    'leave_one_domain_out': loo_domain_beta5,
    'interpretation': 'Drop Fingers (18 rows, 16 resistant) vs Drop Palm (111 rows, 32 resistant) shows how beta5 shifts depending on high-influence clusters.'
}

# 3. Variant-drug records vs unique substitutions and overlap
pos_acv = set(df[(df.gene == 'UL30') & (df.drug == 'Aciclovir')].pos.unique())
pos_fos = set(df[(df.gene == 'UL30') & (df.drug == 'Foscarnet')].pos.unique())
mut_acv = set(zip(df[(df.gene == 'UL30') & (df.drug == 'Aciclovir')].pos, df[(df.gene == 'UL30') & (df.drug == 'Aciclovir')].wt, df[(df.gene == 'UL30') & (df.drug == 'Aciclovir')].mut))
mut_fos = set(zip(df[(df.gene == 'UL30') & (df.drug == 'Foscarnet')].pos, df[(df.gene == 'UL30') & (df.drug == 'Foscarnet')].wt, df[(df.gene == 'UL30') & (df.drug == 'Foscarnet')].mut))

results['variant_record_census'] = {
    'total_variant_drug_records': int(len(df)),
    'total_unique_amino_acid_substitutions': int(len(df.drop_duplicates(subset=['gene', 'pos', 'wt', 'mut']))),
    'ul30_overlap': {
        'acv_positions': len(pos_acv),
        'fos_positions': len(pos_fos),
        'shared_positions': len(pos_acv.intersection(pos_fos)),
        'shared_positions_pct_acv': round(len(pos_acv.intersection(pos_fos)) / len(pos_acv) * 100, 1),
        'acv_substitutions': len(mut_acv),
        'fos_substitutions': len(mut_fos),
        'shared_substitutions': len(mut_acv.intersection(mut_fos)),
        'shared_substitutions_pct_acv': round(len(mut_acv.intersection(mut_fos)) / len(mut_acv) * 100, 1)
    }
}

pathlib.Path('logs').mkdir(exist_ok=True)
with open('logs/step0b_recalibrated_benchmarks.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "="*85)
print("STEP 0b RECALIBRATED BENCHMARKS SUMMARY")
print("="*85)

for k in ['UL30-Aciclovir', 'UL30-Foscarnet', 'UL23-Aciclovir']:
    r = results[k]
    print(f"\nTarget: {k} ({r['description']})")
    print(f"  Evaluable Domains ({r['n_evaluable_domains']}): {r['evaluable_domains_list']}")
    print(f"  Distinct Bootstrap Multisets: {r['distinct_bootstrap_multisets']} (Theoretical Max: {r['theoretical_multisets_max']})")
    print(f"  Domain Means: Pos-Only RF = {r['domain_means']['auc_pos_rf']:.4f} | ESM-2 = {r['domain_means']['auc_esm2']:.4f} | BLOSUM = {r['domain_means']['auc_blosum62']:.4f} | Where = {r['domain_means']['auc_where']:.4f}")
    print(f"  Recalibrated Excess (vs max(0.5, Pos-Only)):")
    print(f"    ESM-2:   {r['excess_metrics']['esm2']['domain_mean_excess']:+.4f} (95% CI: [{r['excess_metrics']['esm2']['ci_95'][0]:+.4f}, {r['excess_metrics']['esm2']['ci_95'][1]:+.4f}]) [Sign: {r['excess_metrics']['esm2']['sign_consistency']}]")
    print(f"    BLOSUM:  {r['excess_metrics']['blosum62']['domain_mean_excess']:+.4f} (95% CI: [{r['excess_metrics']['blosum62']['ci_95'][0]:+.4f}, {r['excess_metrics']['blosum62']['ci_95'][1]:+.4f}]) [Sign: {r['excess_metrics']['blosum62']['sign_consistency']}]")
    print(f"    Where:   {r['excess_metrics']['where']['domain_mean_excess']:+.4f} (95% CI: [{r['excess_metrics']['where']['ci_95'][0]:+.4f}, {r['excess_metrics']['where']['ci_95'][1]:+.4f}]) [Sign: {r['excess_metrics']['where']['sign_consistency']}]")

print("\n" + "="*85)
print("LEAVE-ONE-DOMAIN-OUT BETA5 DIAGNOSTIC (UL30-ACV)")
print("="*85)
for dom, d_res in results['beta5_bimodality_diagnostic']['leave_one_domain_out'].items():
    print(f"  Drop {dom:12s} (N_drop={d_res['n_dropped']:3d}, Pos_drop={d_res['pos_dropped']:2d}): beta5 = {d_res['beta5']:+.4f} (Where beta1 = {d_res['beta1_where']:+.4f})")

print("\n" + "="*85)
print("VARIANT-DRUG RECORD CENSUS & OVERLAP")
print("="*85)
vc = results['variant_record_census']
print(f"  Total Variant-Drug Records: {vc['total_variant_drug_records']}")
print(f"  Total Unique Substitutions: {vc['total_unique_amino_acid_substitutions']}")
print(f"  UL30 ACV vs FOS Shared Positions:     {vc['ul30_overlap']['shared_positions']} / {vc['ul30_overlap']['acv_positions']} ({vc['ul30_overlap']['shared_positions_pct_acv']}%)")
print(f"  UL30 ACV vs FOS Shared Substitutions: {vc['ul30_overlap']['shared_substitutions']} / {vc['ul30_overlap']['acv_substitutions']} ({vc['ul30_overlap']['shared_substitutions_pct_acv']}%)")
