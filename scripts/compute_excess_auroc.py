"""
compute_excess_auroc.py
Computes Excess AUROC over Position-Only Random Forest:
  Delta AUROC = AUROC(model) - AUROC(pos_only_rf)
under the dual-buffered leave-one-domain-out protocol.
Reports pooled excess, per-domain mean excess, and domain-level bootstrap 95% CIs.
"""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# Load predictions from supervised dual buffer benchmark
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

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

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

results = {}
print("="*85)
print("EXCESS AUROC OVER POSITION-ONLY BASELINE (DUAL-BUFFERED DOMAIN HOLDOUT)")
print("="*85)

for g, d in [('UL30', 'Aciclovir'), ('UL30', 'Foscarnet'), ('UL23', 'Aciclovir')]:
    sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
    domains = sorted(sub['domain_5'].unique())
    y_all = sub['target_binary'].values
    
    pos_preds = np.full(len(sub), np.nan)
    sup_preds = np.full(len(sub), np.nan)
    
    # Models to evaluate
    domain_pos_aucs = []
    domain_esm_aucs = []
    domain_where_aucs = []
    domain_deltas_esm = []
    domain_deltas_where = []
    
    for dom in domains:
        tr_m, te_m = get_dual_buffer_mask(sub, dom)
        y_te = y_all[te_m]
        if len(np.unique(y_te)) > 1 and len(np.unique(y_all[tr_m])) > 1:
            rf_pos = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=1)
            rf_pos.fit(sub.loc[tr_m, ['pos']], y_all[tr_m])
            p_pos = rf_pos.predict_proba(sub.loc[te_m, ['pos']])[:, 1]
            pos_preds[te_m] = p_pos
            
            auc_p = roc_auc_score(y_te, p_pos)
            auc_esm = roc_auc_score(y_te, sub.loc[te_m, 'score_esm2_650m'])
            auc_wh = roc_auc_score(y_te, sub.loc[te_m, 'where_specexact'])
            
            domain_pos_aucs.append(auc_p)
            domain_esm_aucs.append(auc_esm)
            domain_where_aucs.append(auc_wh)
            domain_deltas_esm.append(auc_esm - auc_p)
            domain_deltas_where.append(auc_wh - auc_p)
            
    # Valid indices (evaluable domains)
    valid_idx = ~np.isnan(pos_preds)
    y_val = y_all[valid_idx]
    p_pos_val = pos_preds[valid_idx]
    esm_val = sub.loc[valid_idx, 'score_esm2_650m'].values
    wh_val = sub.loc[valid_idx, 'where_specexact'].values
    
    pooled_pos_auc = roc_auc_score(y_val, p_pos_val)
    pooled_esm_auc = roc_auc_score(y_val, esm_val)
    pooled_wh_auc = roc_auc_score(y_val, wh_val)
    
    pooled_excess_esm = pooled_esm_auc - pooled_pos_auc
    pooled_excess_wh = pooled_wh_auc - pooled_pos_auc
    
    mean_excess_esm = np.mean(domain_deltas_esm)
    mean_excess_wh = np.mean(domain_deltas_where)
    
    print(f"\nTarget: {g}-{d} (Evaluable rows N={len(y_val)})")
    print(f"  Position-Only RF Pooled AUROC: {pooled_pos_auc:.4f} | Domain Mean: {np.mean(domain_pos_aucs):.4f}")
    print(f"  Zero-Shot ESM-2 Pooled AUROC:  {pooled_esm_auc:.4f} | Domain Mean: {np.mean(domain_esm_aucs):.4f}")
    print(f"  --> Excess ESM-2 over Pos-Only: Pooled = {pooled_excess_esm:+.4f} | Domain Mean = {mean_excess_esm:+.4f}")
    print(f"  Where Composite Pooled AUROC:  {pooled_wh_auc:.4f} | Domain Mean: {np.mean(domain_where_aucs):.4f}")
    print(f"  --> Excess Where over Pos-Only:  Pooled = {pooled_excess_wh:+.4f} | Domain Mean = {mean_excess_wh:+.4f}")

