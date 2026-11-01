"""
ascertainment_bias_controls.py
Module 0: Ascertainment Bias and Curation Controls for HSV-1 Polymerase (UL30)

Tests whether the high predictive power of 'Where' (AUROC ~0.86 on ACV, ~0.91 on FOS)
reflects genuine functional localization or curation ascertainment bias:
1. Evidence tier stratification: Recombinant / Marker Transfer vs. Clinical Isolates
2. Quantitative assay filtering: Restricting to variants with measured fold-change / IC50
3. Factor decomposition: Evolutionary constraint alone vs. 3D pocket proximity alone
"""
import re
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# 1. Load datasets
df_feat = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
herpes = pd.read_csv('herpesdrg-db/herpesdrg-db.tsv', sep='\t')

# Merge HerpesDRG metadata on unique (gene, aa_change)
herpes_hsv = herpes[herpes['virus'] == 'HSV1'].copy()
# Deduplicate herpesdrg on (gene, aa_change) taking first non-null
herpes_dedup = herpes_hsv.drop_duplicates(subset=['gene', 'aa_change']).copy()

df = pd.merge(df_feat, herpes_dedup[['gene', 'aa_change', 'strain_generation', 'test_method', 'Aciclovir', 'Foscarnet', 'ref_doi']], on=['gene', 'aa_change'], how='left')

# Spec-exact Where
for g in ['UL23', 'UL30']:
    m = (df['gene'] == g)
    z_mu = - (df.loc[m, 'mu_plm_19'] - df.loc[m, 'mu_plm_19'].mean()) / df.loc[m, 'mu_plm_19'].std()
    z_dist = (df.loc[m, 'dist_to_site_angstrom'] - df.loc[m, 'dist_to_site_angstrom'].mean()) / df.loc[m, 'dist_to_site_angstrom'].std()
    df.loc[m, 'where_specexact'] = (z_mu - z_dist) / np.std(z_mu - z_dist)
    df.loc[m, 'constraint_only'] = z_mu
    df.loc[m, 'proximity_only'] = - z_dist # higher = closer

def is_numeric_fc(val):
    if pd.isna(val): return False
    s = str(val).replace('>', '').replace('<', '').strip()
    try:
        float(s)
        return True
    except:
        return False

# Categorize evidence tier
def get_evidence_tier(sg):
    if pd.isna(sg): return 'Unspecified / Literature'
    s = str(sg).strip().lower()
    if 'marker transfer' in s or 'recombinant' in s or 'lab induced' in s or 'lab derived' in s:
        return 'Recombinant / Lab-Validated'
    if 'isolated strain' in s:
        return 'Clinical Isolate'
    return 'Other'

df['evidence_tier'] = [get_evidence_tier(x) for x in df['strain_generation']]

results = {}

print("="*85)
print("MODULE 0: ASCERTAINMENT BIAS & CURATION CONTROL ANALYSIS")
print("="*85)

for drug in ['Aciclovir', 'Foscarnet']:
    sub = df[(df['gene'] == 'UL30') & (df['drug'] == drug)].copy()
    y_full = sub['target_binary'].values
    w_full = sub['where_specexact'].values
    c_full = sub['constraint_only'].values
    p_full = sub['proximity_only'].values
    auc_full_w = roc_auc_score(y_full, w_full)
    auc_full_c = roc_auc_score(y_full, c_full)
    auc_full_p = roc_auc_score(y_full, p_full)
    
    print(f"\nTarget: UL30-{drug} (Full Cohort N={len(sub)}, Pos={y_full.sum()})")
    print(f"  Where Composite AUROC:     {auc_full_w:.4f}")
    print(f"  Constraint Alone AUROC:    {auc_full_c:.4f}")
    print(f"  Pocket Proximity Alone:    {auc_full_p:.4f}")
    
    # 1. Stratified by Evidence Tier
    print(f"\n  --- Control 1: Evidence Tier Stratification ---")
    tier_results = {}
    for tier in ['Clinical Isolate', 'Recombinant / Lab-Validated']:
        t_sub = sub[sub['evidence_tier'] == tier]
        y_t = t_sub['target_binary'].values
        if len(np.unique(y_t)) > 1:
            auc_w = roc_auc_score(y_t, t_sub['where_specexact'].values)
            auc_c = roc_auc_score(y_t, t_sub['constraint_only'].values)
            auc_p = roc_auc_score(y_t, t_sub['proximity_only'].values)
            print(f"  [{tier}] N={len(t_sub)} (Pos={y_t.sum()} / Neg={len(y_t)-y_t.sum()}):")
            print(f"      Where AUROC:        {auc_w:.4f}")
            print(f"      Constraint AUROC:   {auc_c:.4f}")
            print(f"      Proximity AUROC:    {auc_p:.4f}")
            tier_results[tier] = {'N': len(t_sub), 'Pos': int(y_t.sum()), 'auc_where': auc_w, 'auc_constraint': auc_c, 'auc_proximity': auc_p}
        else:
            print(f"  [{tier}] Single-class (N={len(t_sub)}, Pos={y_t.sum()})")
            
    # 2. Control 2: Measured Fold-Change Only Subgroup
    print(f"\n  --- Control 2: Quantitative Fold-Change Measured Subgroup ---")
    sub['has_numeric_fc'] = [is_numeric_fc(v) for v in sub[drug]]
    num_sub = sub[sub['has_numeric_fc']].copy()
    y_num = num_sub['target_binary'].values
    if len(np.unique(y_num)) > 1:
        auc_num_w = roc_auc_score(y_num, num_sub['where_specexact'].values)
        auc_num_c = roc_auc_score(y_num, num_sub['constraint_only'].values)
        auc_num_p = roc_auc_score(y_num, num_sub['proximity_only'].values)
        print(f"  Measured Fold-Change Subset: N={len(num_sub)} (Pos={y_num.sum()} / Neg={len(y_num)-y_num.sum()}):")
        print(f"      Where AUROC:        {auc_num_w:.4f}")
        print(f"      Constraint AUROC:   {auc_num_c:.4f}")
        print(f"      Proximity AUROC:    {auc_num_p:.4f}")
        fc_res = {'N': len(num_sub), 'Pos': int(y_num.sum()), 'auc_where': auc_num_w, 'auc_constraint': auc_num_c, 'auc_proximity': auc_num_p}
    else:
        print(f"  Measured Fold-Change Subset is single-class (N={len(num_sub)})")
        fc_res = None
        
    results[f"UL30-{drug}"] = {
        'full_cohort': {'N': len(sub), 'Pos': int(y_full.sum()), 'auc_where': auc_full_w, 'auc_constraint': auc_full_c, 'auc_proximity': auc_full_p},
        'evidence_tiers': tier_results,
        'measured_fold_change_subset': fc_res
    }

with open('logs/ascertainment_bias_controls.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nSaved ascertainment controls to logs/ascertainment_bias_controls.json")
