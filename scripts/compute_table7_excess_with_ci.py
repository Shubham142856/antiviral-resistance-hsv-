"""
compute_table7_excess_with_ci.py
Computes per-domain AUROC minus Position-Only-RF AUROC:
  Delta_domain = AUROC(model)_d - AUROC(pos_only_rf)_d
Computes domain-clustered bootstrap 95% CIs (1,000 resamples of evaluable domains).
Also computes position-clustered bootstrap 95% CIs for pooled excess.
Saves results to logs/table7_excess_auroc_ci.json.
"""
import json
import pathlib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier

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

targets = [('UL30', 'Aciclovir'), ('UL30', 'Foscarnet'), ('UL23', 'Aciclovir')]
results = {}

print("="*90)
print("TABLE 7: EXCESS AUROC OVER POSITION-ONLY RF WITH BOOTSTRAP CONFIDENCE INTERVALS")
print("="*90)

rng = np.random.RandomState(42)

for g, d in targets:
    sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
    domains = sorted(sub['domain_5'].unique())
    y_all = sub['target_binary'].values
    
    pos_preds = np.full(len(sub), np.nan)
    
    domain_records = []
    
    for dom in domains:
        tr_m, te_m = get_dual_buffer_mask(sub, dom)
        y_te = y_all[te_m]
        # Check evaluability: at least one of each class in train and test
        if len(np.unique(y_te)) > 1 and len(np.unique(y_all[tr_m])) > 1:
            rf_pos = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=1)
            rf_pos.fit(sub.loc[tr_m, ['pos']], y_all[tr_m])
            p_pos = rf_pos.predict_proba(sub.loc[te_m, ['pos']])[:, 1]
            pos_preds[te_m] = p_pos
            
            auc_p = roc_auc_score(y_te, p_pos)
            auc_esm = roc_auc_score(y_te, sub.loc[te_m, 'score_esm2_650m'])
            auc_bert = roc_auc_score(y_te, sub.loc[te_m, 'score_protbert'])
            auc_blo = roc_auc_score(y_te, sub.loc[te_m, 'score_blosum62'])
            auc_wh = roc_auc_score(y_te, sub.loc[te_m, 'where_specexact'])
            
            domain_records.append({
                'domain': dom,
                'n_test': len(y_te),
                'n_pos': int(y_te.sum()),
                'auc_pos': auc_p,
                'auc_esm': auc_esm,
                'auc_bert': auc_bert,
                'auc_blo': auc_blo,
                'auc_wh': auc_wh,
                'delta_esm': auc_esm - auc_p,
                'delta_bert': auc_bert - auc_p,
                'delta_blo': auc_blo - auc_p,
                'delta_wh': auc_wh - auc_p
            })
            
    # Domain-clustered bootstrap over evaluable domains
    n_dom = len(domain_records)
    boot_domain_esm = []
    boot_domain_wh = []
    boot_domain_bert = []
    boot_domain_blo = []
    
    for _ in range(2000):
        resampled_doms = rng.choice(domain_records, size=n_dom, replace=True)
        boot_domain_esm.append(np.mean([r['delta_esm'] for r in resampled_doms]))
        boot_domain_wh.append(np.mean([r['delta_wh'] for r in resampled_doms]))
        boot_domain_bert.append(np.mean([r['delta_bert'] for r in resampled_doms]))
        boot_domain_blo.append(np.mean([r['delta_blo'] for r in resampled_doms]))
        
    ci_esm_dom = [float(np.percentile(boot_domain_esm, 2.5)), float(np.percentile(boot_domain_esm, 97.5))]
    ci_wh_dom = [float(np.percentile(boot_domain_wh, 2.5)), float(np.percentile(boot_domain_wh, 97.5))]
    ci_bert_dom = [float(np.percentile(boot_domain_bert, 2.5)), float(np.percentile(boot_domain_bert, 97.5))]
    ci_blo_dom = [float(np.percentile(boot_domain_blo, 2.5)), float(np.percentile(boot_domain_blo, 97.5))]
    
    # Valid indices for pooled
    valid_idx = ~np.isnan(pos_preds)
    y_val = y_all[valid_idx]
    p_pos_val = pos_preds[valid_idx]
    esm_val = sub.loc[valid_idx, 'score_esm2_650m'].values
    wh_val = sub.loc[valid_idx, 'where_specexact'].values
    bert_val = sub.loc[valid_idx, 'score_protbert'].values
    blo_val = sub.loc[valid_idx, 'score_blosum62'].values
    
    pooled_pos = roc_auc_score(y_val, p_pos_val)
    pooled_esm = roc_auc_score(y_val, esm_val)
    pooled_wh = roc_auc_score(y_val, wh_val)
    pooled_bert = roc_auc_score(y_val, bert_val)
    pooled_blo = roc_auc_score(y_val, blo_val)
    
    # Position-clustered bootstrap on evaluable subset for pooled excess
    val_sub = sub.loc[valid_idx].copy().reset_index(drop=True)
    u_pos = val_sub['pos'].unique()
    pos_map = {p: val_sub[val_sub['pos'] == p].index.values for p in u_pos}
    
    boot_pooled_esm = []
    boot_pooled_wh = []
    boot_pooled_bert = []
    boot_pooled_blo = []
    
    for _ in range(1000):
        samp_p = rng.choice(u_pos, size=len(u_pos), replace=True)
        samp_idx = np.concatenate([pos_map[p] for p in samp_p])
        y_b = y_val[samp_idx]
        if len(np.unique(y_b)) > 1:
            auc_pos_b = roc_auc_score(y_b, p_pos_val[samp_idx])
            boot_pooled_esm.append(roc_auc_score(y_b, esm_val[samp_idx]) - auc_pos_b)
            boot_pooled_wh.append(roc_auc_score(y_b, wh_val[samp_idx]) - auc_pos_b)
            boot_pooled_bert.append(roc_auc_score(y_b, bert_val[samp_idx]) - auc_pos_b)
            boot_pooled_blo.append(roc_auc_score(y_b, blo_val[samp_idx]) - auc_pos_b)
            
    ci_pooled_esm = [float(np.percentile(boot_pooled_esm, 2.5)), float(np.percentile(boot_pooled_esm, 97.5))]
    ci_pooled_wh = [float(np.percentile(boot_pooled_wh, 2.5)), float(np.percentile(boot_pooled_wh, 97.5))]
    ci_pooled_bert = [float(np.percentile(boot_pooled_bert, 2.5)), float(np.percentile(boot_pooled_bert, 97.5))]
    ci_pooled_blo = [float(np.percentile(boot_pooled_blo, 2.5)), float(np.percentile(boot_pooled_blo, 97.5))]
    
    mean_excess_esm = float(np.mean([r['delta_esm'] for r in domain_records]))
    mean_excess_wh = float(np.mean([r['delta_wh'] for r in domain_records]))
    mean_excess_bert = float(np.mean([r['delta_bert'] for r in domain_records]))
    mean_excess_blo = float(np.mean([r['delta_blo'] for r in domain_records]))
    
    split_valid = (0.45 <= pooled_pos <= 0.55)
    
    results[f"{g}-{d}"] = {
        'n_evaluable_rows': int(len(y_val)),
        'n_evaluable_domains': n_dom,
        'domain_records': domain_records,
        'pos_only_rf_pooled': float(pooled_pos),
        'split_validity_confirmed': bool(split_valid),
        'excess_esm2': {
            'domain_mean_excess': mean_excess_esm,
            'domain_bootstrap_95ci': ci_esm_dom,
            'pooled_excess': float(pooled_esm - pooled_pos),
            'pooled_position_bootstrap_95ci': ci_pooled_esm
        },
        'excess_where': {
            'domain_mean_excess': mean_excess_wh,
            'domain_bootstrap_95ci': ci_wh_dom,
            'pooled_excess': float(pooled_wh - pooled_pos),
            'pooled_position_bootstrap_95ci': ci_pooled_wh
        },
        'excess_protbert': {
            'domain_mean_excess': mean_excess_bert,
            'domain_bootstrap_95ci': ci_bert_dom,
            'pooled_excess': float(pooled_bert - pooled_pos),
            'pooled_position_bootstrap_95ci': ci_pooled_bert
        },
        'excess_blosum62': {
            'domain_mean_excess': mean_excess_blo,
            'domain_bootstrap_95ci': ci_blo_dom,
            'pooled_excess': float(pooled_blo - pooled_pos),
            'pooled_position_bootstrap_95ci': ci_pooled_blo
        }
    }
    
    print(f"\nTarget: {g}-{d} (Evaluable Domains: {n_dom}, Evaluable Rows: {len(y_val)})")
    print(f"  Pos-Only RF Pooled AUROC: {pooled_pos:.4f} -> Split Validity in [0.45, 0.55]: {split_valid}")
    print(f"  ESM-2 Excess AUROC:")
    print(f"    Per-Domain Mean Excess: {mean_excess_esm:+.4f} (Domain 95% CI: [{ci_esm_dom[0]:+.4f}, {ci_esm_dom[1]:+.4f}])")
    print(f"    Pooled Excess:          {pooled_esm - pooled_pos:+.4f} (Pos-Cluster 95% CI: [{ci_pooled_esm[0]:+.4f}, {ci_pooled_esm[1]:+.4f}])")
    print(f"  Where Composite Excess AUROC:")
    print(f"    Per-Domain Mean Excess: {mean_excess_wh:+.4f} (Domain 95% CI: [{ci_wh_dom[0]:+.4f}, {ci_wh_dom[1]:+.4f}])")
    print(f"    Pooled Excess:          {pooled_wh - pooled_pos:+.4f} (Pos-Cluster 95% CI: [{ci_pooled_wh[0]:+.4f}, {ci_pooled_wh[1]:+.4f}])")

pathlib.Path('logs').mkdir(exist_ok=True)
with open('logs/table7_excess_auroc_ci.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nSaved Table 7 excess AUROC results to logs/table7_excess_auroc_ci.json")

