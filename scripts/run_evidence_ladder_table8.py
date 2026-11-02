"""
run_evidence_ladder_table8.py
=============================
Computes Table 8: The Evidence Ladder across 4 designs:
  - D1: Pooled full cohort (out-of-fold CV for Position-RF / pooled zero-shot with 95% bootstrap CI)
  - D2: Dual-buffer domain-mean (holding out structural domains with 12A + 5aa buffer)
  - D3: S0 within-study pair-weighted (unified study-clustered 95% bootstrap CI with seed 42) & study-equal mean
  - S2: Numeric fold-change matched-assay stratum (balanced quantitative phenotypes)

Targets:
  - UL30-Aciclovir (Primary)
  - UL23-Aciclovir (Primary)
  - UL30-Foscarnet (Descriptive)

Models evaluated:
  - Position-Only RF (D1 out-of-fold and D2 out-of-domain)
  - Where (Positional Prior: spec-exact Z(Constraint) - Z(Dist))
  - ESM-2 650M
  - ProtBERT
  - -entropy (Conservation)
  - PSSM
  - BLOSUM62
"""

import sys
import json
import warnings
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
import re

warnings.filterwarnings('ignore')

print("Loading data for Evidence Ladder...", flush=True)
df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
spatial_df = pd.read_csv('data/processed/primary_3d_spatial_folds.tsv', sep='\t')
cmi_all = pd.read_csv('data/processed/cmi_2025_all_sheets.csv')
herpes = pd.read_csv('herpesdrg-db/herpesdrg-db.tsv', sep='\t')

cols_to_merge = [
    'coord_x', 'coord_y', 'coord_z', 'has_3d_coord',
    'score_esm2_650m', 'score_protbert', 'score_neg_entropy',
    'score_pssm', 'score_blosum62'
]
df = pd.merge(
    df,
    spatial_df[['gene', 'pos', 'wt', 'mut', 'drug'] + cols_to_merge],
    on=['gene', 'pos', 'wt', 'mut', 'drug'],
    how='left'
)

# Compute spec-exact Where
for g in ['UL23', 'UL30']:
    m = (df['gene'] == g)
    z_mu = - (df.loc[m, 'mu_plm_19'] - df.loc[m, 'mu_plm_19'].mean()) / df.loc[m, 'mu_plm_19'].std()
    z_dist = (df.loc[m, 'dist_to_site_angstrom'] - df.loc[m, 'dist_to_site_angstrom'].mean()) / df.loc[m, 'dist_to_site_angstrom'].std()
    df.loc[m, 'where_specexact'] = (z_mu - z_dist) / np.std(z_mu - z_dist)

def parse_mut(s):
    match = re.match(r'^([A-Z])(\d+)([A-Z])$', str(s).strip())
    if match:
        return match.group(1), int(match.group(2)), match.group(3)
    return None, None, None

hsv1_herpes = herpes[herpes['virus'] == 'HSV1'].copy()
parsed = [parse_mut(x) for x in hsv1_herpes['aa_change']]
hsv1_herpes['wt'] = [p[0] for p in parsed]
hsv1_herpes['pos'] = [p[1] for p in parsed]
hsv1_herpes['mut'] = [p[2] for p in parsed]
hsv1_herpes['gene_std'] = hsv1_herpes['gene'].astype(str).str.strip().replace({'TK': 'UL23', 'POL': 'UL30', 'DPOL': 'UL30'})
hsv1_herpes_dedup = hsv1_herpes.drop_duplicates(subset=['gene_std', 'pos', 'wt', 'mut']).copy()

def get_hdrg_drug_type(row):
    val = row['Aciclovir'] if row['drug'] == 'Aciclovir' else row['Foscarnet']
    if pd.isna(val) or str(val).strip() in ['nan', 'None', '']:
        return '<Missing in HerpesDRG>'
    s = str(val).replace('>', '').replace('<', '').strip()
    try:
        float(s)
        return 'Numeric Fold-Change (e.g. 2.5, 5.0)'
    except:
        return f"Categorical String: '{val.strip()}'"

df = pd.merge(
    df,
    hsv1_herpes_dedup[['gene_std', 'pos', 'wt', 'mut', 'Aciclovir', 'Foscarnet', 'ref_doi', 'ref_title']],
    left_on=['gene', 'pos', 'wt', 'mut'],
    right_on=['gene_std', 'pos', 'wt', 'mut'],
    how='left'
)
df['herpesdrg_phenotype_format'] = [get_hdrg_drug_type(r) for _, r in df.iterrows()]

cmi_hsv1 = cmi_all[cmi_all['Source_Sheet'].isin(['HSV-1 TK', 'HSV-1 Pol'])].copy()
cmi_hsv1['gene_std'] = cmi_hsv1['Gene'].astype(str).str.strip().replace({'Thymidine Kinase': 'UL23', 'DNA Polymerase': 'UL30'})
c_parsed = [parse_mut(x) for x in cmi_hsv1['AA change']]
cmi_hsv1['wt'] = [p[0] for p in c_parsed]
cmi_hsv1['pos'] = [p[1] for p in c_parsed]
cmi_hsv1['mut'] = [p[2] for p in c_parsed]
cmi_hsv1['drug_std'] = cmi_hsv1['Drug'].astype(str).str.strip()
cmi_hsv1_dedup = cmi_hsv1.drop_duplicates(subset=['gene_std', 'drug_std', 'pos', 'wt', 'mut']).copy()

df = pd.merge(
    df,
    cmi_hsv1_dedup[['gene_std', 'drug_std', 'pos', 'wt', 'mut', 'PMIDs']],
    left_on=['gene', 'drug', 'pos', 'wt', 'mut'],
    right_on=['gene_std', 'drug_std', 'pos', 'wt', 'mut'],
    how='left'
)

def get_study_id(row):
    doi = str(row['ref_doi']).strip()
    if doi and doi not in ['nan', 'None', '']:
        return f"DOI:{doi}"
    pmid = str(row['PMIDs']).strip()
    if pmid and pmid not in ['nan', 'None', '']:
        return f"PMID:{pmid}"
    title = str(row['ref_title']).strip()
    if title and title not in ['nan', 'None', '']:
        return f"Title:{title[:40]}"
    return "Unassigned"

df['study_id'] = [get_study_id(r) for _, r in df.iterrows()]
print(f"Data merged: N={len(df)}", flush=True)

# Fast AUC implementation
def fast_auc(y_true, y_score):
    n1 = int(np.sum(y_true == 1))
    n0 = len(y_true) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    ranks = rankdata(y_score)
    u = np.sum(ranks[y_true == 1]) - (n1 * (n1 + 1)) / 2.0
    return float(u / (n1 * n0))

def fast_bootstrap_ci(y_true, y_score, n_boot=2000, seed=42):
    rng = np.random.RandomState(seed)
    n = len(y_true)
    boot_aucs = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, size=n)
        yt_b = y_true[idx]
        n1 = int(np.sum(yt_b == 1))
        n0 = len(yt_b) - n1
        if n1 == 0 or n0 == 0:
            continue
        ranks = rankdata(y_score[idx])
        u = np.sum(ranks[yt_b == 1]) - (n1 * (n1 + 1)) / 2.0
        boot_aucs.append(u / (n1 * n0))
    if len(boot_aucs) == 0:
        return np.nan, np.nan
    return float(np.percentile(boot_aucs, 2.5)), float(np.percentile(boot_aucs, 97.5))

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

SCORE_NAMES = [
    ('Position-Only RF (Out-of-Fold)', 'pos_rf'),
    ('Where (Positional Prior)', 'where_specexact'),
    ('ESM-2 650M', 'score_esm2_650m'),
    ('ProtBERT', 'score_protbert'),
    ('-entropy', 'score_neg_entropy'),
    ('PSSM', 'score_pssm'),
    ('BLOSUM62', 'score_blosum62')
]

TARGETS = [
    ('UL30', 'Aciclovir', 'Primary Target', ['Exo', 'Palm', 'Fingers', 'Thumb_Cterm']),
    ('UL30', 'Foscarnet', 'Descriptive / Secondary', ['Palm', 'Fingers']),
    ('UL23', 'Aciclovir', 'Primary Target (Comparative Kinase)', ['ATP_Core', 'Nucleoside_Pocket', 'Lid_Dimer', 'C-term'])
]

# -------------------------------------------------------------
# D3: Unified Study-Clustered Bootstrap setup (matching Step 1)
# -------------------------------------------------------------
print("Setting up unified study-clustered bootstrap for D3 (S0)...", flush=True)

# Build per-study pairs and concordances for each target and each score
target_s0_data = {}
all_mixed_studies = set()

for g, d, desc, _ in TARGETS:
    t_key = f"{g}-{d}"
    sub = df[(df['gene'] == g) & (df['drug'] == d)].copy()
    st_groups = sub.groupby('study_id')
    
    t_scores_data = {}
    for s_label, s_col in SCORE_NAMES:
        if s_col == 'pos_rf':
            continue
        s_stats = {}
        for st_id, grp in st_groups:
            rv = grp.loc[grp['target_binary'] == 1, s_col].values.astype(float)
            sv = grp.loc[grp['target_binary'] == 0, s_col].values.astype(float)
            nr = len(rv)
            ns = len(sv)
            if nr > 0 and ns > 0:
                n_p = nr * ns
                diffs = rv[:, None] - sv[None, :]
                conc = float(np.sum(diffs > 0) + 0.5 * np.sum(diffs == 0))
                s_stats[st_id] = {'pairs': n_p, 'conc': conc, 'auc': conc / n_p}
                all_mixed_studies.add(st_id)
        t_scores_data[s_col] = s_stats
    target_s0_data[t_key] = t_scores_data

study_universe = sorted(list(all_mixed_studies))
n_studies = len(study_universe)
study_to_idx = {st: i for i, st in enumerate(study_universe)}
print(f"Total mixed studies in unified universe: {n_studies}", flush=True)

# Precompute pair and conc arrays of shape (n_studies,)
pair_arrays = {}
conc_arrays = {s_col: {} for _, s_col in SCORE_NAMES if s_col != 'pos_rf'}

for g, d, desc, _ in TARGETS:
    t_key = f"{g}-{d}"
    p_arr = np.zeros(n_studies)
    for st_id, s_data in target_s0_data[t_key]['where_specexact'].items():
        p_arr[study_to_idx[st_id]] = s_data['pairs']
    pair_arrays[t_key] = p_arr
    for _, s_col in SCORE_NAMES:
        if s_col == 'pos_rf':
            continue
        c_arr = np.zeros(n_studies)
        for st_id, s_data in target_s0_data[t_key][s_col].items():
            c_arr[study_to_idx[st_id]] = s_data['conc']
        conc_arrays[s_col][t_key] = c_arr

# Run unified study-clustered bootstrap with seed=42 (identical to Step 1)
np.random.seed(42)
B = 2000
boot_d3 = {s_col: {f"{g}-{d}": [] for g, d, _, _ in TARGETS} for _, s_col in SCORE_NAMES if s_col != 'pos_rf'}

for b in range(B):
    resampled_idx = np.random.choice(n_studies, size=n_studies, replace=True)
    counts = np.bincount(resampled_idx, minlength=n_studies)
    for g, d, _, _ in TARGETS:
        t_key = f"{g}-{d}"
        p_tot = np.dot(counts, pair_arrays[t_key])
        for _, s_col in SCORE_NAMES:
            if s_col == 'pos_rf':
                continue
            if p_tot > 0:
                auc_b = float(np.dot(counts, conc_arrays[s_col][t_key]) / p_tot)
                boot_d3[s_col][t_key].append(auc_b)
            else:
                boot_d3[s_col][t_key].append(np.nan)

# Compute point estimates and CIs for D3
d3_final = {}
for g, d, desc, _ in TARGETS:
    t_key = f"{g}-{d}"
    d3_final[t_key] = {}
    for s_label, s_col in SCORE_NAMES:
        if s_col == 'pos_rf':
            d3_final[t_key][s_col] = {'pair_weighted': np.nan, 'study_mean': np.nan, 'ci_95': (np.nan, np.nan)}
            continue
        tot_pairs = np.sum(pair_arrays[t_key])
        tot_conc = np.sum(conc_arrays[s_col][t_key])
        pw_auc = float(tot_conc / tot_pairs) if tot_pairs > 0 else np.nan
        st_aucs = [v['auc'] for v in target_s0_data[t_key][s_col].values()]
        sm_auc = float(np.mean(st_aucs)) if len(st_aucs) > 0 else np.nan
        
        clean = [x for x in boot_d3[s_col][t_key] if not np.isnan(x)]
        ci_pw = (float(np.percentile(clean, 2.5)), float(np.percentile(clean, 97.5))) if len(clean) > 0 else (np.nan, np.nan)
        d3_final[t_key][s_col] = {
            'pair_weighted': pw_auc,
            'study_mean': sm_auc,
            'ci_95': ci_pw
        }

table8_results = {}

for g, d, desc, evaluable_domains in TARGETS:
    t_key = f"{g}-{d}"
    print(f"\nProcessing {t_key} ({desc})...", flush=True)
    sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
    y = sub['target_binary'].values.astype(int)
    n_tot = len(sub)
    n_res = int(y.sum())
    n_sens = n_tot - n_res
    
    t_res = {'N_total': n_tot, 'N_res': n_res, 'N_sens': n_sens, 'scores': {}}
    
    # ---------------------------------------------------------
    # D1: Pooled Full Cohort
    # ---------------------------------------------------------
    print(f"  Computing D1 (Pooled)...", flush=True)
    # Position RF (out-of-fold)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    rf_oof = np.zeros(n_tot)
    for tr_idx, te_idx in skf.split(sub, y):
        rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5)
        rf.fit(sub.iloc[tr_idx][['pos']], y[tr_idx])
        rf_oof[te_idx] = rf.predict_proba(sub.iloc[te_idx][['pos']])[:, 1]
    rf_d1_auc = fast_auc(y, rf_oof)
    rf_d1_ci = fast_bootstrap_ci(y, rf_oof, seed=42)
    
    d1_dict = {'pos_rf': {'auc': rf_d1_auc, 'ci_95': rf_d1_ci}}
    for s_label, s_col in SCORE_NAMES:
        if s_col == 'pos_rf':
            continue
        vals = sub[s_col].values.astype(float)
        auc_v = fast_auc(y, vals)
        ci_v = fast_bootstrap_ci(y, vals, seed=42)
        d1_dict[s_col] = {'auc': auc_v, 'ci_95': ci_v}
        
    # ---------------------------------------------------------
    # D2: Dual-Buffer Domain-Mean
    # ---------------------------------------------------------
    print(f"  Computing D2 (Dual-Buffer Domain-Mean)...", flush=True)
    d2_per_domain = {s_col: [] for _, s_col in SCORE_NAMES}
    for dom in evaluable_domains:
        tr_m, te_m = get_dual_buffer_mask(sub, dom)
        y_te = y[te_m]
        sub_te = sub[te_m]
        
        # Position RF trained out-of-domain
        rf_d2 = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5)
        rf_d2.fit(sub.loc[tr_m, ['pos']], y[tr_m])
        te_preds = rf_d2.predict_proba(sub.loc[te_m, ['pos']])[:, 1]
        if len(np.unique(te_preds)) < 2:
            rf_d2_auc = 0.5000
        else:
            rf_d2_auc = fast_auc(y_te, te_preds)
        d2_per_domain['pos_rf'].append(rf_d2_auc)
        
        for s_label, s_col in SCORE_NAMES:
            if s_col == 'pos_rf':
                continue
            vals_te = sub_te[s_col].values.astype(float)
            d2_per_domain[s_col].append(fast_auc(y_te, vals_te))
            
    d2_dict = {}
    rng_d = np.random.RandomState(42)
    n_doms = len(evaluable_domains)
    for s_label, s_col in SCORE_NAMES:
        arr = np.array(d2_per_domain[s_col])
        m_auc = float(np.mean(arr))
        if n_doms >= 4:
            b_means = [float(np.mean(rng_d.choice(arr, size=n_doms, replace=True))) for _ in range(2000)]
            ci_d = (float(np.percentile(b_means, 2.5)), float(np.percentile(b_means, 97.5)))
        else:
            ci_d = None
        d2_dict[s_col] = {'auc': m_auc, 'ci_95': ci_d, 'per_domain': [float(x) for x in arr]}
        
    # ---------------------------------------------------------
    # D3: S0 Within-Study Pairing (From Unified Bootstrap)
    # ---------------------------------------------------------
    d3_dict = d3_final[t_key]
        
    # ---------------------------------------------------------
    # S2: Numeric Fold-Change Stratum
    # ---------------------------------------------------------
    print(f"  Computing S2 (Numeric FC Stratum)...", flush=True)
    s2_sub = sub[sub['herpesdrg_phenotype_format'] == 'Numeric Fold-Change (e.g. 2.5, 5.0)'].copy().reset_index(drop=True)
    y_s2 = s2_sub['target_binary'].values.astype(int)
    n_s2 = len(s2_sub)
    n_s2_r = int(y_s2.sum())
    n_s2_s = n_s2 - n_s2_r
    print(f"    S2 N={n_s2} (Res={n_s2_r}, Sens={n_s2_s})", flush=True)
    
    s2_dict = {}
    for s_label, s_col in SCORE_NAMES:
        if s_col == 'pos_rf':
            s2_dict[s_col] = {'auc': np.nan, 'ci_95': (np.nan, np.nan)}
            continue
        vals_s2 = s2_sub[s_col].values.astype(float)
        auc_s2 = fast_auc(y_s2, vals_s2)
        ci_s2 = fast_bootstrap_ci(y_s2, vals_s2, seed=42)
        s2_dict[s_col] = {'auc': auc_s2, 'ci_95': ci_s2}
        
    for s_label, s_col in SCORE_NAMES:
        t_res['scores'][s_label] = {
            'D1_pooled': d1_dict[s_col],
            'D2_domain_mean': d2_dict[s_col],
            'D3_within_study': d3_dict[s_col],
            'S2_numeric_fc': s2_dict[s_col]
        }
    table8_results[t_key] = t_res

with open('logs/evidence_ladder_table8_results.json', 'w') as f:
    json.dump(table8_results, f, indent=2)

print("\nSUCCESS: All calculations completed and saved to logs/evidence_ladder_table8_results.json", flush=True)
