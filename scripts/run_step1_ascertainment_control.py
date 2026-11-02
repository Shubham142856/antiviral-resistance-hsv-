"""
run_step1_ascertainment_control.py
==================================
Step 1: Within-Study Ascertainment Control and Source Diagnostics.
Anchored strictly to frozen pre-registration: docs/STEP1_ASCERTAINMENT_CONTROL_PREREGISTRATION.md
(Commit: ef69018, SHA-256: 61500c501df95c28bf44692b91b805d3f1eb8bbcdc37b0fcbc6abe0d7ea5b435)

Computes:
  A. Positivity Check:
     - Mixed-study count (M)
     - Total same-study R-S pairs
     - Top 5 studies by pair count (absolute and percentage share)
     - Distinct resistant variants contributing to same-study pairs
     - Positivity rules check (>=100 pairs, >=5 mixed studies, <=50% single study share)
  B. S0 (Primary):
     - Pair-weighted AUROC on same-study pairs
     - Sensitivity: Study-equal-weighted mean of per-study AUROCs
     - Leave-one-study-out (LOSO) S0 AUROC for Where and BLOSUM
     - Study-clustered bootstrap (B=2000) with shared study resample for cross-gene contrast D
  C. Secondary Strata (Descriptive only, no pass/fail):
     - S1: CMI cell-culture confirmed = yes
     - S2: HerpesDRG numeric fold-change
     - S3: Marker transfer + BAC (UL30-ACV, flagged underpowered)
  D. Source Diagnostics (Descriptive, within label class):
     - Within sensitive rows: Where and BLOSUM AUROC for survey-only vs mixed studies; in-CMI vs not
     - Within resistant rows: same if n>=10 in both bins
  E. Pre-Registered Decision Verdicts (Where survival and Dissociation survival)

Outputs:
  - logs/step1_ascertainment_control_results.json
"""

import os
import re
import json
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')

# 1. Load Primary Benchmark Cohort
primary_df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')

# 2. Spec-Exact Where and BLOSUM Features
for g in ['UL23', 'UL30']:
    m = (primary_df['gene'] == g)
    z_mu = - (primary_df.loc[m, 'mu_plm_19'] - primary_df.loc[m, 'mu_plm_19'].mean()) / primary_df.loc[m, 'mu_plm_19'].std()
    z_dist = (primary_df.loc[m, 'dist_to_site_angstrom'] - primary_df.loc[m, 'dist_to_site_angstrom'].mean()) / primary_df.loc[m, 'dist_to_site_angstrom'].std()
    primary_df.loc[m, 'where_specexact'] = (z_mu - z_dist) / np.std(z_mu - z_dist)

# Verify BLOSUM column
primary_df['blosum_score'] = primary_df['severity_blosum'].values

# 3. Load Literature Metadata (HerpesDRG and CMI)
herpes = pd.read_csv('herpesdrg-db/herpesdrg-db.tsv', sep='\t')
hsv1_herpes = herpes[herpes['virus'] == 'HSV1'].copy()

def parse_mut(s):
    match = re.match(r'^([A-Z])(\d+)([A-Z])$', str(s).strip())
    if match:
        return match.group(1), int(match.group(2)), match.group(3)
    return None, None, None

parsed = [parse_mut(x) for x in hsv1_herpes['aa_change']]
hsv1_herpes['wt'] = [p[0] for p in parsed]
hsv1_herpes['pos'] = [p[1] for p in parsed]
hsv1_herpes['mut'] = [p[2] for p in parsed]
hsv1_herpes['gene_std'] = hsv1_herpes['gene'].astype(str).str.strip().replace({'TK': 'UL23', 'POL': 'UL30', 'DPOL': 'UL30'})

for col in ['strain_generation', 'test_method', 'control_strain']:
    hsv1_herpes[col] = hsv1_herpes[col].astype(str).str.strip()
    hsv1_herpes.loc[hsv1_herpes[col].isin(['nan', 'None', '']), col] = np.nan

hdrg_cols = ['strain_generation', 'test_method', 'control_strain', 'Aciclovir', 'Foscarnet', 'ref_doi', 'ref_title']
hsv1_herpes_dedup = hsv1_herpes.drop_duplicates(subset=['gene_std', 'pos', 'wt', 'mut']).copy()

# Load CMI
cmi = pd.read_csv('data/processed/cmi_2025_all_sheets.csv')
cmi_hsv1 = cmi[cmi['Source_Sheet'].isin(['HSV-1 TK', 'HSV-1 Pol'])].copy()
cmi_hsv1['gene_std'] = cmi_hsv1['Gene'].astype(str).str.strip().replace({'Thymidine Kinase': 'UL23', 'DNA Polymerase': 'UL30'})
cmi_parsed = [parse_mut(x) for x in cmi_hsv1['AA change']]
cmi_hsv1['wt'] = [p[0] for p in cmi_parsed]
cmi_hsv1['pos'] = [p[1] for p in cmi_parsed]
cmi_hsv1['mut'] = [p[2] for p in cmi_parsed]
cmi_hsv1['drug_std'] = cmi_hsv1['Drug'].astype(str).str.strip()

for col in ['resistance status confirmed in cell culture', 'clinical susceptibility information available', 'All genes sequenced that may confer resistance?', 'Resistance']:
    cmi_hsv1[col] = cmi_hsv1[col].astype(str).str.strip()
    cmi_hsv1.loc[cmi_hsv1[col].isin(['nan', 'None', '']), col] = np.nan

cmi_cols = [
    'resistance status confirmed in cell culture',
    'clinical susceptibility information available',
    'All genes sequenced that may confer resistance?',
    'Resistance',
    'PMIDs'
]
cmi_hsv1_dedup = cmi_hsv1.drop_duplicates(subset=['gene_std', 'drug_std', 'pos', 'wt', 'mut']).copy()

# Merge onto primary_df
merged = pd.merge(
    primary_df,
    hsv1_herpes_dedup[['gene_std', 'pos', 'wt', 'mut'] + hdrg_cols],
    left_on=['gene', 'pos', 'wt', 'mut'],
    right_on=['gene_std', 'pos', 'wt', 'mut'],
    how='left'
)

merged = pd.merge(
    merged,
    cmi_hsv1_dedup[['gene_std', 'drug_std', 'pos', 'wt', 'mut'] + cmi_cols],
    left_on=['gene', 'drug', 'pos', 'wt', 'mut'],
    right_on=['gene_std', 'drug_std', 'pos', 'wt', 'mut'],
    how='left'
)

# Helper to check if drug value in HerpesDRG is numeric
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

merged['herpesdrg_phenotype_format'] = [get_hdrg_drug_type(r) for _, r in merged.iterrows()]

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

merged['study_id'] = [get_study_id(r) for _, r in merged.iterrows()]

# Target definitions
targets = [
    ('UL30', 'Aciclovir', 'Primary'),
    ('UL30', 'Foscarnet', 'Secondary'),
    ('UL23', 'Aciclovir', 'Counter-Arm')
]

# Helper function for Mann-Whitney U AUROC on pairs
def get_study_pair_stats(df_sub, score_col):
    """
    Computes per-study concordant pairs, total pairs, and AUROC.
    """
    st_groups = df_sub.groupby('study_id')
    stats = {}
    for st_id, grp in st_groups:
        res_scores = grp.loc[grp['target_binary'] == 1, score_col].values
        sens_scores = grp.loc[grp['target_binary'] == 0, score_col].values
        n_r = len(res_scores)
        n_s = len(sens_scores)
        if n_r > 0 and n_s > 0:
            n_p = n_r * n_s
            diffs = res_scores[:, None] - sens_scores[None, :]
            concordant = float(np.sum(diffs > 0) + 0.5 * np.sum(diffs == 0))
            stats[st_id] = {
                'n_res': n_r,
                'n_sens': n_s,
                'pairs': n_p,
                'concordant': concordant,
                'auc': concordant / n_p
            }
    return stats

def aggregate_study_stats(stats_dict):
    total_pairs = sum(s['pairs'] for s in stats_dict.values())
    total_concordant = sum(s['concordant'] for s in stats_dict.values())
    pw_auc = (total_concordant / total_pairs) if total_pairs > 0 else np.nan
    study_equal_auc = float(np.mean([s['auc'] for s in stats_dict.values()])) if len(stats_dict) > 0 else np.nan
    return pw_auc, study_equal_auc, total_pairs

print("="*95)
print("STEP 1: WITHIN-STUDY ASCERTAINMENT CONTROL & SOURCE DIAGNOSTICS")
print("Anchored to pre-registration: docs/STEP1_ASCERTAINMENT_CONTROL_PREREGISTRATION.md")
print("="*95)

# A. POSITIVITY CHECK & S0 COMPUTATION
positivity_summary = {}
study_pair_data = {}  # {target: {'where': stats, 'blosum': stats}}

for g, d, role in targets:
    t_name = f"{g}-{d}"
    sub = merged[(merged['gene'] == g) & (merged['drug'] == d)].copy().reset_index(drop=True)
    
    st_summary = sub.groupby('study_id')['target_binary'].agg(
        total='count',
        res=lambda x: int((x == 1).sum()),
        sens=lambda x: int((x == 0).sum())
    )
    mixed_st = st_summary[(st_summary['res'] > 0) & (st_summary['sens'] > 0)].copy()
    mixed_st['pairs'] = mixed_st['res'] * mixed_st['sens']
    mixed_st = mixed_st.sort_values(by='pairs', ascending=False)
    
    total_pairs = int(mixed_st['pairs'].sum())
    mixed_count = len(mixed_st)
    
    top5 = []
    for st_id, row in mixed_st.head(5).iterrows():
        p_cnt = int(row['pairs'])
        share = float(p_cnt / total_pairs * 100) if total_pairs > 0 else 0.0
        top5.append({
            'study_id': st_id,
            'resistant': int(row['res']),
            'sensitive': int(row['sens']),
            'pairs': p_cnt,
            'share_pct': round(share, 2)
        })
        
    top1_share = top5[0]['share_pct'] if top5 else 0.0
    
    # Distinct resistant variants in mixed studies (Amendment Item 1)
    mixed_studies_set = set(mixed_st.index)
    res_in_mixed = sub[sub['study_id'].isin(mixed_studies_set) & (sub['target_binary'] == 1)]
    distinct_res_vars = res_in_mixed[['wt', 'pos', 'mut']].drop_duplicates()
    distinct_res_count = len(distinct_res_vars)
    
    rule_pairs = (total_pairs >= 100)
    rule_mixed = (mixed_count >= 5)
    rule_share = (top1_share <= 50.0)
    passes_positivity = rule_pairs and rule_mixed and rule_share
    
    pos_status = "PASS" if passes_positivity else "INCONCLUSIVE (Thin Data / Concentration)"
    if t_name == "UL30-Foscarnet":
        pos_status = "THIN DATA (Secondary Only, Non-independent)"
        
    positivity_summary[t_name] = {
        'role': role,
        'mixed_studies_count': mixed_count,
        'total_same_study_pairs': total_pairs,
        'top1_study_share_pct': top1_share,
        'distinct_resistant_variants': distinct_res_count,
        'total_resistant_in_mixed': len(res_in_mixed),
        'top5_studies': top5,
        'rule_pairs_ge_100': rule_pairs,
        'rule_mixed_ge_5': rule_mixed,
        'rule_share_le_50': rule_share,
        'positivity_verdict': pos_status
    }
    
    # Compute precalculated per-study stats
    study_pair_data[t_name] = {
        'where': get_study_pair_stats(sub, 'where_specexact'),
        'blosum': get_study_pair_stats(sub, 'blosum_score')
    }

# B. S0 EVALUATION: POINT ESTIMATES AND LEAVE-ONE-STUDY-OUT
s0_results = {}

for g, d, role in targets:
    t_name = f"{g}-{d}"
    st_w = study_pair_data[t_name]['where']
    st_b = study_pair_data[t_name]['blosum']
    
    pw_where, me_where, n_pairs = aggregate_study_stats(st_w)
    pw_blosum, me_blosum, _ = aggregate_study_stats(st_b)
    
    # Leave-One-Study-Out S0 AUROC (Amendment Item 1)
    loso_where = {}
    loso_blosum = {}
    for omitted_study in st_w.keys():
        sub_w = {k: v for k, v in st_w.items() if k != omitted_study}
        sub_b = {k: v for k, v in st_b.items() if k != omitted_study}
        pw_w_loso, _, _ = aggregate_study_stats(sub_w)
        pw_b_loso, _, _ = aggregate_study_stats(sub_b)
        loso_where[omitted_study] = pw_w_loso
        loso_blosum[omitted_study] = pw_b_loso
        
    s0_results[t_name] = {
        'role': role,
        'pairs': n_pairs,
        'where': {
            'pair_weighted_auroc': pw_where,
            'study_equal_auroc': me_where,
            'loso_auroc': loso_where
        },
        'blosum62': {
            'pair_weighted_auroc': pw_blosum,
            'study_equal_auroc': me_blosum,
            'loso_auroc': loso_blosum
        },
        'delta_where_minus_blosum': (pw_where - pw_blosum) if not np.isnan(pw_where) and not np.isnan(pw_blosum) else np.nan
    }

# C. FAST VECTORIZED STUDY-CLUSTERED BOOTSTRAP (B=2000) WITH SHARED RESAMPLING FOR CONTRAST D
print("\nRunning fast vectorized study-clustered bootstrap (B=2000 resamples)...")
np.random.seed(42)
B = 2000

# Unified universe of mixed studies across all targets
all_mixed_studies = set()
for t_name in s0_results:
    all_mixed_studies.update(study_pair_data[t_name]['where'].keys())
study_universe = sorted(list(all_mixed_studies))
n_studies = len(study_universe)
study_to_idx = {st: i for i, st in enumerate(study_universe)}

# Precompute arrays for each target: shape (n_studies,)
pair_arrays = {}
conc_w_arrays = {}
conc_b_arrays = {}

for t_name in s0_results:
    p_arr = np.zeros(n_studies)
    cw_arr = np.zeros(n_studies)
    cb_arr = np.zeros(n_studies)
    for st_id, s_data in study_pair_data[t_name]['where'].items():
        idx = study_to_idx[st_id]
        p_arr[idx] = s_data['pairs']
        cw_arr[idx] = s_data['concordant']
    for st_id, s_data in study_pair_data[t_name]['blosum'].items():
        idx = study_to_idx[st_id]
        cb_arr[idx] = s_data['concordant']
    pair_arrays[t_name] = p_arr
    conc_w_arrays[t_name] = cw_arr
    conc_b_arrays[t_name] = cb_arr

boot_where = {t_name: [] for t_name in s0_results}
boot_blosum = {t_name: [] for t_name in s0_results}
boot_diff = {t_name: [] for t_name in s0_results}
boot_contrast_d = []

for b in range(B):
    # Resample study indices with replacement
    resampled_idx = np.random.choice(n_studies, size=n_studies, replace=True)
    counts = np.bincount(resampled_idx, minlength=n_studies)
    
    b_vals = {}
    for t_name in s0_results:
        p_tot = np.dot(counts, pair_arrays[t_name])
        if p_tot > 0:
            w_auc = float(np.dot(counts, conc_w_arrays[t_name]) / p_tot)
            b_auc = float(np.dot(counts, conc_b_arrays[t_name]) / p_tot)
            boot_where[t_name].append(w_auc)
            boot_blosum[t_name].append(b_auc)
            diff = w_auc - b_auc
            boot_diff[t_name].append(diff)
            b_vals[t_name] = diff
        else:
            boot_where[t_name].append(np.nan)
            boot_blosum[t_name].append(np.nan)
            boot_diff[t_name].append(np.nan)
            
    if 'UL30-Aciclovir' in b_vals and 'UL23-Aciclovir' in b_vals:
        boot_contrast_d.append(b_vals['UL30-Aciclovir'] - b_vals['UL23-Aciclovir'])
    else:
        boot_contrast_d.append(np.nan)

# Compute 95% Percentile Bootstrap Intervals
for t_name in s0_results:
    w_clean = [x for x in boot_where[t_name] if not np.isnan(x)]
    b_clean = [x for x in boot_blosum[t_name] if not np.isnan(x)]
    d_clean = [x for x in boot_diff[t_name] if not np.isnan(x)]
    
    ci_w = [float(np.percentile(w_clean, 2.5)), float(np.percentile(w_clean, 97.5))] if len(w_clean) > 0 else [np.nan, np.nan]
    ci_b = [float(np.percentile(b_clean, 2.5)), float(np.percentile(b_clean, 97.5))] if len(b_clean) > 0 else [np.nan, np.nan]
    ci_d = [float(np.percentile(d_clean, 2.5)), float(np.percentile(d_clean, 97.5))] if len(d_clean) > 0 else [np.nan, np.nan]
    
    s0_results[t_name]['where']['bootstrap_95ci'] = ci_w
    s0_results[t_name]['blosum62']['bootstrap_95ci'] = ci_b
    s0_results[t_name]['delta_where_minus_blosum_95ci'] = ci_d

d_clean_all = [x for x in boot_contrast_d if not np.isnan(x)]
contrast_d_point = s0_results['UL30-Aciclovir']['delta_where_minus_blosum'] - s0_results['UL23-Aciclovir']['delta_where_minus_blosum']
contrast_d_ci = [float(np.percentile(d_clean_all, 2.5)), float(np.percentile(d_clean_all, 97.5))] if len(d_clean_all) > 0 else [np.nan, np.nan]

contrast_d_results = {
    'point_estimate': contrast_d_point,
    'bootstrap_95ci': contrast_d_ci,
    'n_resamples': len(d_clean_all)
}

# D. SECONDARY STRATA (DESCRIPTIVE ONLY, NO PASS/FAIL)
secondary_strata_results = {}

for g, d, role in targets:
    t_name = f"{g}-{d}"
    sub = merged[(merged['gene'] == g) & (merged['drug'] == d)].copy().reset_index(drop=True)
    t_dict = {}
    
    # S1: CMI cell culture confirmed == 'yes'
    s1_mask = (sub['resistance status confirmed in cell culture'] == 'yes')
    s1_sub = sub[s1_mask]
    s1_r = int(s1_sub['target_binary'].sum())
    s1_s = len(s1_sub) - s1_r
    if s1_r > 0 and s1_s > 0:
        auc_w = float(roc_auc_score(s1_sub['target_binary'], s1_sub['where_specexact']))
        auc_b = float(roc_auc_score(s1_sub['target_binary'], s1_sub['blosum_score']))
    else:
        auc_w, auc_b = np.nan, np.nan
    t_dict['S1_CMI_Confirmed'] = {
        'N': len(s1_sub),
        'Resistant': s1_r,
        'Sensitive': s1_s,
        'evaluable': min(s1_r, s1_s) >= 10,
        'where_auroc': auc_w,
        'blosum_auroc': auc_b
    }
    
    # S2: HerpesDRG numeric fold-change
    s2_mask = (sub['herpesdrg_phenotype_format'] == 'Numeric Fold-Change (e.g. 2.5, 5.0)')
    s2_sub = sub[s2_mask]
    s2_r = int(s2_sub['target_binary'].sum())
    s2_s = len(s2_sub) - s2_r
    if s2_r > 0 and s2_s > 0:
        auc_w = float(roc_auc_score(s2_sub['target_binary'], s2_sub['where_specexact']))
        auc_b = float(roc_auc_score(s2_sub['target_binary'], s2_sub['blosum_score']))
    else:
        auc_w, auc_b = np.nan, np.nan
    t_dict['S2_Numeric_FC'] = {
        'N': len(s2_sub),
        'Resistant': s2_r,
        'Sensitive': s2_s,
        'evaluable': min(s2_r, s2_s) >= 10,
        'where_auroc': auc_w,
        'blosum_auroc': auc_b
    }
    
    # S3: Marker transfer + BAC (UL30-ACV descriptive row)
    if t_name == 'UL30-Aciclovir':
        s3_mask = sub['strain_generation'].isin(['Marker Transfer', 'Recombinant BAC', 'Marker transfer'])
        s3_sub = sub[s3_mask]
        s3_r = int(s3_sub['target_binary'].sum())
        s3_s = len(s3_sub) - s3_r
        if s3_r > 0 and s3_s > 0:
            auc_w = float(roc_auc_score(s3_sub['target_binary'], s3_sub['where_specexact']))
            auc_b = float(roc_auc_score(s3_sub['target_binary'], s3_sub['blosum_score']))
        else:
            auc_w, auc_b = np.nan, np.nan
        t_dict['S3_Recombinant_Transfer_BAC'] = {
            'N': len(s3_sub),
            'Resistant': s3_r,
            'Sensitive': s3_s,
            'evaluable': False,  # Flagged underpowered as pre-registered
            'flag': 'Underpowered (N=37, Res=22, Sens=15)',
            'where_auroc': auc_w,
            'blosum_auroc': auc_b
        }
        
    secondary_strata_results[t_name] = t_dict

# E. SOURCE DIAGNOSTICS (WITHIN LABEL CLASS, DESCRIPTIVE ONLY)
source_diagnostics = {}

for g, d, role in targets:
    t_name = f"{g}-{d}"
    sub = merged[(merged['gene'] == g) & (merged['drug'] == d)].copy().reset_index(drop=True)
    
    st_counts = sub.groupby('study_id')['target_binary'].agg(
        res=lambda x: int((x == 1).sum()),
        sens=lambda x: int((x == 0).sum())
    )
    mixed_st_set = set(st_counts[(st_counts['res'] > 0) & (st_counts['sens'] > 0)].index)
    
    sub['is_mixed_study'] = sub['study_id'].isin(mixed_st_set).astype(int)
    sub['in_cmi'] = sub['resistance status confirmed in cell culture'].notna().astype(int)
    
    # 1. Within Sensitive Rows (y=0)
    sens_sub = sub[sub['target_binary'] == 0].copy()
    diag_sens = {}
    
    vc_mix = sens_sub['is_mixed_study'].value_counts().to_dict()
    if len(vc_mix) > 1 and min(vc_mix.values()) >= 5:
        auc_w_mix = float(roc_auc_score(sens_sub['is_mixed_study'], sens_sub['where_specexact']))
        auc_b_mix = float(roc_auc_score(sens_sub['is_mixed_study'], sens_sub['blosum_score']))
        diag_sens['survey_only_vs_mixed'] = {
            'counts': vc_mix,
            'evaluable': True,
            'where_auroc': auc_w_mix,
            'blosum_auroc': auc_b_mix
        }
    else:
        diag_sens['survey_only_vs_mixed'] = {'counts': vc_mix, 'evaluable': False}
        
    vc_cmi = sens_sub['in_cmi'].value_counts().to_dict()
    if len(vc_cmi) > 1 and min(vc_cmi.values()) >= 5:
        auc_w_cmi = float(roc_auc_score(sens_sub['in_cmi'], sens_sub['where_specexact']))
        auc_b_cmi = float(roc_auc_score(sens_sub['in_cmi'], sens_sub['blosum_score']))
        diag_sens['in_cmi_vs_not'] = {
            'counts': vc_cmi,
            'evaluable': True,
            'where_auroc': auc_w_cmi,
            'blosum_auroc': auc_b_cmi
        }
    else:
        diag_sens['in_cmi_vs_not'] = {'counts': vc_cmi, 'evaluable': False}
        
    # 2. Within Resistant Rows (y=1)
    res_sub = sub[sub['target_binary'] == 1].copy()
    diag_res = {}
    
    vc_res_mix = res_sub['is_mixed_study'].value_counts().to_dict()
    if len(vc_res_mix) > 1 and min(vc_res_mix.values()) >= 10:
        auc_w_rmix = float(roc_auc_score(res_sub['is_mixed_study'], res_sub['where_specexact']))
        auc_b_rmix = float(roc_auc_score(res_sub['is_mixed_study'], res_sub['blosum_score']))
        diag_res['mixed_vs_all_resistant'] = {
            'counts': vc_res_mix,
            'evaluable': True,
            'where_auroc': auc_w_rmix,
            'blosum_auroc': auc_b_rmix
        }
    else:
        diag_res['mixed_vs_all_resistant'] = {
            'counts': vc_res_mix,
            'evaluable': False,
            'note': f"Underpowered (counts={vc_res_mix}, min < 10)"
        }
        
    vc_res_cmi = res_sub['in_cmi'].value_counts().to_dict()
    if len(vc_res_cmi) > 1 and min(vc_res_cmi.values()) >= 10:
        auc_w_rcmi = float(roc_auc_score(res_sub['in_cmi'], res_sub['where_specexact']))
        auc_b_rcmi = float(roc_auc_score(res_sub['in_cmi'], res_sub['blosum_score']))
        diag_res['in_cmi_vs_not'] = {
            'counts': vc_res_cmi,
            'evaluable': True,
            'where_auroc': auc_w_rcmi,
            'blosum_auroc': auc_b_rcmi
        }
    else:
        diag_res['in_cmi_vs_not'] = {
            'counts': vc_res_cmi,
            'evaluable': False,
            'note': f"Underpowered (counts={vc_res_cmi}, min < 10)"
        }
        
    source_diagnostics[t_name] = {
        'within_sensitive': diag_sens,
        'within_resistant': diag_res
    }

# F. PRE-REGISTERED DECISION VERDICTS
# 1. UL30-ACV Where Survival Test:
#    PASS: S0 pair-weighted point >= 0.65 AND CI lower > 0.55
#    ARTIFACT: point < 0.55
#    INCONCLUSIVE: otherwise
ul30_pw = s0_results['UL30-Aciclovir']['where']['pair_weighted_auroc']
ul30_ci = s0_results['UL30-Aciclovir']['where']['bootstrap_95ci']

if ul30_pw >= 0.65 and ul30_ci[0] > 0.55:
    where_verdict = "PASS"
elif ul30_pw < 0.55:
    where_verdict = "ARTIFACT"
else:
    where_verdict = "INCONCLUSIVE"

# 2. Dissociation Survival Test:
#    PASS: D > 0 and CI lower > 0
#    INCONCLUSIVE: otherwise
d_point = contrast_d_results['point_estimate']
d_ci = contrast_d_results['bootstrap_95ci']

if d_point > 0 and d_ci[0] > 0:
    dissociation_verdict = "PASS"
else:
    dissociation_verdict = "INCONCLUSIVE"

verdicts = {
    'UL30_ACV_Where_Survival': {
        'point_estimate': ul30_pw,
        'ci_95': ul30_ci,
        'thresholds': {'pass_point_ge': 0.65, 'pass_ci_lower_gt': 0.55, 'artifact_point_lt': 0.55},
        'verdict': where_verdict
    },
    'Cross_Gene_Dissociation_Contrast_D': {
        'point_estimate': d_point,
        'ci_95': d_ci,
        'rule': 'PASS if D > 0 and 95% CI excludes 0',
        'verdict': dissociation_verdict
    }
}

# Compile complete results object
final_output = {
    'protocol_metadata': {
        'protocol_version': '1.0 (Audit-Grade Frozen Specification with Pre-Run Amendment)',
        'preregistration_doc': 'docs/STEP1_ASCERTAINMENT_CONTROL_PREREGISTRATION.md',
        'amended_commit_hash': 'ef69018dfa2ea2d545783233a8518cb71f7d870d',
        'amended_sha256': '61500c501df95c28bf44692b91b805d3f1eb8bbcdc37b0fcbc6abe0d7ea5b435'
    },
    'positivity_check': positivity_summary,
    'S0_primary_within_study': s0_results,
    'cross_gene_contrast_D': contrast_d_results,
    'secondary_strata_descriptive': secondary_strata_results,
    'source_diagnostics_within_class': source_diagnostics,
    'decision_verdicts': verdicts
}

# Save to logs
os.makedirs('logs', exist_ok=True)
with open('logs/step1_ascertainment_control_results.json', 'w') as f:
    json.dump(final_output, f, indent=2)

print("\n" + "="*95)
print("EXECUTION COMPLETED SUCCESSFULLY. RESULTS SAVED TO logs/step1_ascertainment_control_results.json")
print("="*95)
