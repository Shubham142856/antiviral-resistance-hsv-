"""
compute_step0b_final.py
Step 0b Final Benchmark & Audit:
1. Table 7: Apply pre-registered N_pos >= 5 rule.
   - UL30-ACV: Exo, Palm, Fingers, Thumb (n=4)
   - UL30-FOS: Palm, Fingers (n=2, dropping Exo [N_pos=2] and Thumb [N_pos=3])
   - UL23-ACV: ATP Core, Nucleoside Pocket, Lid/Dimer, C-term (n=4)
2. Per-domain table: N, N_pos, AUROC for each model:
   - Position-Only RF
   - Supervised LR (trained out-of-domain with L2 penalty, features: ESM-2, ProtBERT, ESM2-35M, Cons, BLOSUM62, Where, Sev)
   - ESM-2 650M
   - ProtBERT
   - BLOSUM62
   - Where (spec-exact: Z(Constraint) - Z(Dist))
   - Excess = AUROC - max(0.50, AUROC_pos_only)
3. Domain-bootstrap:
   - Print number of distinct domain-resamples per target.
   - For targets with < 4 evaluable domains (FOS, n=2), report NO bootstrap CI, only per-domain values and mean.
   - For 4-domain targets (UL30-ACV, UL23-ACV), CIs are descriptive only.
4. Excess reference = max(0.50, position-only RF domain mean). Recompute all excess values and BLOSUM62 CIs.
   - Show raw domain means against 0.50 chance line.
   - Report UL23 minus UL30 BLOSUM62 domain-mean difference as descriptive.
5. Verify UL23 BLOSUM62 domain mean comes from domain folds and print per-domain values.
6. Primary claim for Where = domain-mean values; pooled reported as secondary.
7. Saves results to logs/step0b_final_recalibrated_benchmarks.json.
"""
import json
import math
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings('ignore')

# 1. Load data
df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
spatial_df = pd.read_csv('data/processed/primary_3d_spatial_folds.tsv', sep='\t')
cols_to_merge = ['coord_x', 'coord_y', 'coord_z', 'has_3d_coord', 'score_esm2_650m', 'score_protbert', 'score_blosum62', 'score_weighted_cons', 'score_esm2_35m']
df = pd.merge(df, spatial_df[['gene', 'pos', 'wt', 'mut', 'drug'] + cols_to_merge], on=['gene', 'pos', 'wt', 'mut', 'drug'], how='left')

# Spec-exact Where and Severity
for g in ['UL23', 'UL30']:
    m = (df['gene'] == g)
    z_mu = - (df.loc[m, 'mu_plm_19'] - df.loc[m, 'mu_plm_19'].mean()) / df.loc[m, 'mu_plm_19'].std()
    z_dist = (df.loc[m, 'dist_to_site_angstrom'] - df.loc[m, 'dist_to_site_angstrom'].mean()) / df.loc[m, 'dist_to_site_angstrom'].std()
    df.loc[m, 'where_specexact'] = (z_mu - z_dist) / np.std(z_mu - z_dist)
    df.loc[m, 'sev_specexact'] = - (df.loc[m, 'severity_plm'] - df.loc[m, 'severity_plm'].mean()) / df.loc[m, 'severity_plm'].std()

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
    ('UL30', 'Aciclovir', 'Primary Target', ['Exo', 'Palm', 'Fingers', 'Thumb_Cterm']),
    ('UL30', 'Foscarnet', 'Secondary Target (Strict N_pos>=5 Rule: Palm, Fingers only)', ['Palm', 'Fingers']),
    ('UL23', 'Aciclovir', 'Comparative Kinase (Dispensable Target)', ['ATP_Core', 'Nucleoside_Pocket', 'Lid_Dimer', 'C-term'])
]

sup_features = ['score_esm2_650m', 'score_protbert', 'score_esm2_35m', 'score_weighted_cons', 'score_blosum62', 'where_specexact', 'sev_specexact']

results = {}
rng = np.random.RandomState(42)

print("="*115)
print("FINAL STEP 0b BENCHMARK: SUPERVISED LR, ZERO-SHOT MODELS & RECALIBRATED EXCESS METRICS")
print("="*115)

for g, d, desc, evaluable_domains in targets:
    sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
    y_all = sub['target_binary'].values
    
    # Fill NAs in features
    for col in sup_features:
        if sub[col].isna().any():
            sub[col] = sub[col].fillna(sub[col].median())
            
    domain_records = []
    
    print(f"\nTarget: {g}-{d} | {desc}")
    print("-" * 115)
    print(f"{'Domain':<18} | {'N (Res)':<9} | {'Pos-RF':<7} | {'RefBaseline':<11} | {'ESM-2':<7} | {'ProtBERT':<8} | {'BLOSUM62':<8} | {'Where':<7} | {'Sup-LR':<7}")
    print("-" * 115)
    
    for dom in evaluable_domains:
        tr_m, te_m = get_dual_buffer_mask(sub, dom)
        y_te = y_all[te_m]
        n_pos_test = int(y_te.sum())
        
        # Position-only RF
        rf_pos = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=1)
        rf_pos.fit(sub.loc[tr_m, ['pos']], y_all[tr_m])
        p_pos = rf_pos.predict_proba(sub.loc[te_m, ['pos']])[:, 1]
        auc_p = float(roc_auc_score(y_te, p_pos))
        
        # Zero-shot models
        auc_esm = float(roc_auc_score(y_te, sub.loc[te_m, 'score_esm2_650m']))
        auc_bert = float(roc_auc_score(y_te, sub.loc[te_m, 'score_protbert']))
        auc_blo = float(roc_auc_score(y_te, sub.loc[te_m, 'score_blosum62']))
        auc_wh = float(roc_auc_score(y_te, sub.loc[te_m, 'where_specexact']))
        
        # Supervised Logistic Regression L2
        X_tr = sub.loc[tr_m, sup_features].values
        X_te = sub.loc[te_m, sup_features].values
        mu_tr = np.mean(X_tr, axis=0)
        std_tr = np.std(X_tr, axis=0)
        std_tr[std_tr == 0] = 1.0
        X_tr_s = (X_tr - mu_tr) / std_tr
        X_te_s = (X_te - mu_tr) / std_tr
        
        lr = LogisticRegression(C=1.0, max_iter=1000, random_state=42)
        lr.fit(X_tr_s, y_all[tr_m])
        p_lr = lr.predict_proba(X_te_s)[:, 1]
        auc_lr = float(roc_auc_score(y_te, p_lr))
        
        # Reference baseline: max(0.50, auc_pos_rf)
        ref_baseline = max(0.50, auc_p)
        
        rec = {
            'domain': dom,
            'n_test': int(len(y_te)),
            'n_pos': n_pos_test,
            'auc_pos_rf': auc_p,
            'reference_baseline': ref_baseline,
            'auc_esm2': auc_esm,
            'auc_protbert': auc_bert,
            'auc_blosum62': auc_blo,
            'auc_where': auc_wh,
            'auc_sup_lr': auc_lr,
            'excess_esm2': auc_esm - ref_baseline,
            'excess_protbert': auc_bert - ref_baseline,
            'excess_blosum62': auc_blo - ref_baseline,
            'excess_where': auc_wh - ref_baseline,
            'excess_sup_lr': auc_lr - ref_baseline,
            'raw_diff_05_esm2': auc_esm - 0.50,
            'raw_diff_05_blosum62': auc_blo - 0.50,
            'raw_diff_05_where': auc_wh - 0.50,
            'raw_diff_05_sup_lr': auc_lr - 0.50
        }
        domain_records.append(rec)
        print(f"{dom:<18} | {rec['n_test']} ({rec['n_pos']})    | {auc_p:.4f}  | {ref_baseline:.4f}      | {auc_esm:.4f}  | {auc_bert:.4f}   | {auc_blo:.4f}   | {auc_wh:.4f}  | {auc_lr:.4f}")
        
    n_dom = len(domain_records)
    
    # Calculate domain means
    mean_pos = float(np.mean([r['auc_pos_rf'] for r in domain_records]))
    ref_mean = max(0.50, mean_pos)
    mean_esm = float(np.mean([r['auc_esm2'] for r in domain_records]))
    mean_bert = float(np.mean([r['auc_protbert'] for r in domain_records]))
    mean_blo = float(np.mean([r['auc_blosum62'] for r in domain_records]))
    mean_wh = float(np.mean([r['auc_where'] for r in domain_records]))
    mean_lr = float(np.mean([r['auc_sup_lr'] for r in domain_records]))
    
    excess_mean_esm = mean_esm - ref_mean
    excess_mean_bert = mean_bert - ref_mean
    excess_mean_blo = mean_blo - ref_mean
    excess_mean_wh = mean_wh - ref_mean
    excess_mean_lr = mean_lr - ref_mean
    
    # Sign consistency vs position-only RF
    sign_esm = f"{sum(r['auc_esm2'] > r['auc_pos_rf'] for r in domain_records)}/{n_dom}"
    sign_bert = f"{sum(r['auc_protbert'] > r['auc_pos_rf'] for r in domain_records)}/{n_dom}"
    sign_blo = f"{sum(r['auc_blosum62'] > r['auc_pos_rf'] for r in domain_records)}/{n_dom}"
    sign_wh = f"{sum(r['auc_where'] > r['auc_pos_rf'] for r in domain_records)}/{n_dom}"
    sign_lr = f"{sum(r['auc_sup_lr'] > r['auc_pos_rf'] for r in domain_records)}/{n_dom}"
    
    # Distinct multiset count: comb(n_dom + n_dom - 1, n_dom)
    theoretical_distinct = math.comb(2 * n_dom - 1, n_dom)
    
    distinct_multisets = set()
    boot_excess_esm = []
    boot_excess_wh = []
    boot_excess_blo = []
    boot_excess_lr = []
    
    if n_dom >= 4:
        for _ in range(2000):
            idx_resamp = tuple(sorted(rng.choice(range(n_dom), size=n_dom, replace=True)))
            distinct_multisets.add(idx_resamp)
            boot_excess_esm.append(np.mean([domain_records[i]['excess_esm2'] for i in idx_resamp]))
            boot_excess_wh.append(np.mean([domain_records[i]['excess_where'] for i in idx_resamp]))
            boot_excess_blo.append(np.mean([domain_records[i]['excess_blosum62'] for i in idx_resamp]))
            boot_excess_lr.append(np.mean([domain_records[i]['excess_sup_lr'] for i in idx_resamp]))
            
        ci_esm = [float(np.percentile(boot_excess_esm, 2.5)), float(np.percentile(boot_excess_esm, 97.5))]
        ci_wh = [float(np.percentile(boot_excess_wh, 2.5)), float(np.percentile(boot_excess_wh, 97.5))]
        ci_blo = [float(np.percentile(boot_excess_blo, 2.5)), float(np.percentile(boot_excess_blo, 97.5))]
        ci_lr = [float(np.percentile(boot_excess_lr, 2.5)), float(np.percentile(boot_excess_lr, 97.5))]
        ci_note = "Descriptive only (small-n domain bootstrap)"
    else:
        # Strict rule: NO bootstrap CI for targets with < 4 evaluable domains
        ci_esm = None
        ci_wh = None
        ci_blo = None
        ci_lr = None
        ci_note = f"NO bootstrap CI reported (< 4 evaluable domains, n={n_dom}, theoretical distinct resamples = {theoretical_distinct})"
        
    print("-" * 115)
    print(f"{'Domain Mean':<18} | {'Total':<9} | {mean_pos:.4f}  | {ref_mean:.4f}      | {mean_esm:.4f}  | {mean_bert:.4f}   | {mean_blo:.4f}   | {mean_wh:.4f}  | {mean_lr:.4f}")
    print(f"{'Excess (vs Ref)':<18} | {'--':<9} | {'0.0000':<7} | {'--':<11} | {excess_mean_esm:+.4f}  | {excess_mean_bert:+.4f}   | {excess_mean_blo:+.4f}   | {excess_mean_wh:+.4f}  | {excess_mean_lr:+.4f}")
    print(f"{'Sign Consistency':<18} | {'--':<9} | {'--':<7} | {'--':<11} | {sign_esm:<7}  | {sign_bert:<8}   | {sign_blo:<8}   | {sign_wh:<7}  | {sign_lr:<7}")
    print(f"Distinct Resamples: {theoretical_distinct} | CI Status: {ci_note}")
    if ci_wh is not None:
        print(f"  Where Excess 95% CI: [{ci_wh[0]:.4f}, {ci_wh[1]:.4f}]")
        print(f"  ESM-2 Excess 95% CI: [{ci_esm[0]:.4f}, {ci_esm[1]:.4f}]")
        print(f"  BLOSUM62 Excess 95% CI: [{ci_blo[0]:.4f}, {ci_blo[1]:.4f}]")
        print(f"  Sup-LR Excess 95% CI: [{ci_lr[0]:.4f}, {ci_lr[1]:.4f}]")
        
    results[f"{g}-{d}"] = {
        'description': desc,
        'n_evaluable_domains': n_dom,
        'evaluable_domains_list': evaluable_domains,
        'theoretical_distinct_resamples': theoretical_distinct,
        'ci_status': ci_note,
        'per_domain_records': domain_records,
        'domain_means': {
            'pos_only_rf': mean_pos,
            'reference_clamped': ref_mean,
            'esm2_650m': mean_esm,
            'protbert': mean_bert,
            'blosum62': mean_blo,
            'where_specexact': mean_wh,
            'supervised_lr': mean_lr
        },
        'excess_over_reference': {
            'esm2_650m': excess_mean_esm,
            'protbert': excess_mean_bert,
            'blosum62': excess_mean_blo,
            'where_specexact': excess_mean_wh,
            'supervised_lr': excess_mean_lr
        },
        'excess_descriptive_ci_95': {
            'esm2_650m': ci_esm,
            'where_specexact': ci_wh,
            'blosum62': ci_blo,
            'supervised_lr': ci_lr
        },
        'sign_consistency_vs_pos_rf': {
            'esm2_650m': sign_esm,
            'protbert': sign_bert,
            'blosum62': sign_blo,
            'where_specexact': sign_wh,
            'supervised_lr': sign_lr
        }
    }

# Compute descriptive difference in BLOSUM62 domain-mean: UL23 minus UL30-ACV
blo_ul23_mean = results['UL23-Aciclovir']['domain_means']['blosum62']
blo_ul30_mean = results['UL30-Aciclovir']['domain_means']['blosum62']
delta_blo = blo_ul23_mean - blo_ul30_mean

print("\n" + "="*115)
print("DESCRIPTIVE DISSOCIATION CONTRAST: BLOSUM62 RADICALNESS ACROSS GENES")
print("="*115)
print(f"UL23-Aciclovir BLOSUM62 Domain Mean: {blo_ul23_mean:.4f} (Excess over 0.50: {blo_ul23_mean - 0.50:+.4f})")
print(f"UL30-Aciclovir BLOSUM62 Domain Mean: {blo_ul30_mean:.4f} (Excess over 0.50: {blo_ul30_mean - 0.50:+.4f})")
print(f"Descriptive Difference (UL23 - UL30): Delta = {delta_blo:+.4f}")
print("Label: DESCRIPTIVE ONLY (evaluated across n=4 structural domains per gene; no inferential claim)")

results['descriptive_contrast_blosum62'] = {
    'ul23_mean': blo_ul23_mean,
    'ul30_mean': blo_ul30_mean,
    'delta': delta_blo,
    'label': 'DESCRIPTIVE ONLY (n=4 domains each)'
}

with open('logs/step0b_final_recalibrated_benchmarks.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nSaved Final Step 0b benchmark to logs/step0b_final_recalibrated_benchmarks.json")
