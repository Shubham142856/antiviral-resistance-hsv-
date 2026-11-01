"""
evaluate_supervised_dual_buffer.py
Evaluates supervised models (Logistic Regression, HistGradientBoosting)
vs Zero-shot models (ESM-2, ProtBERT, BLOSUM62) and Position-Only RF
under the validated dual-buffered Leave-One-Domain-Out split (seq_buf=5, spatial_buf=12.0 A).
Anchored to Phase 2 Item 5 of the evaluation protocol.
"""

import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

print("Loading datasets...", flush=True)
df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
spatial_df = pd.read_csv('data/processed/primary_3d_spatial_folds.tsv', sep='\t')

# Merge 3D coordinates and scores
cols_to_merge = ['coord_x', 'coord_y', 'coord_z', 'has_3d_coord', 'score_blosum62', 'score_neg_entropy', 'score_weighted_cons', 'score_pssm', 'score_protbert', 'score_esm2_35m', 'score_esm2_650m']
df = pd.merge(df, spatial_df[['gene', 'pos', 'wt', 'mut', 'drug'] + cols_to_merge], on=['gene', 'pos', 'wt', 'mut', 'drug'], how='left')

# Pre-reg signed features
df['neg_esm2_650m'] = df['score_esm2_650m'] # higher = deleterious
df['neg_protbert'] = df['score_protbert']
df['neg_esm2_35m'] = df['score_esm2_35m']
df['weighted_conservation'] = df['score_weighted_cons']
df['blosum_radicalness'] = df['score_blosum62']

# Pre-reg signed Where: Z(-mu_plm_19) - Z(dist)
for g in ['UL23', 'UL30']:
    m = (df['gene'] == g)
    z_mu = - (df.loc[m, 'mu_plm_19'] - df.loc[m, 'mu_plm_19'].mean()) / df.loc[m, 'mu_plm_19'].std()
    z_dist = (df.loc[m, 'dist_to_site_angstrom'] - df.loc[m, 'dist_to_site_angstrom'].mean()) / df.loc[m, 'dist_to_site_angstrom'].std()
    df.loc[m, 'where_prereg'] = (z_mu - z_dist) / np.std(z_mu - z_dist)
    df.loc[m, 'sev_prereg'] = - (df.loc[m, 'severity_plm'] - df.loc[m, 'severity_plm'].mean()) / df.loc[m, 'severity_plm'].std()

def create_dual_buffer_split(sub_df, test_domain, seq_buf=5, spatial_buf=12.0):
    te_mask = (sub_df['domain_5'] == test_domain)
    te_idx = sub_df[te_mask].index.values
    test_positions = set(sub_df.loc[te_idx, 'pos'].unique())
    
    # Collect test 3D coords
    has_coord_test = sub_df[te_mask & (sub_df['has_3d_coord'] == 1)]
    if len(has_coord_test) > 0:
        test_cas = has_coord_test[['coord_x', 'coord_y', 'coord_z']].drop_duplicates().values
    else:
        test_cas = np.empty((0, 3))
    
    cand_tr = sub_df[~te_mask]
    drop_tr_idx = []
    
    for idx_val, r in cand_tr.iterrows():
        p = int(r['pos'])
        # 1. Sequence buffer exclusion
        if min([abs(p - tp) for tp in test_positions]) <= seq_buf:
            drop_tr_idx.append(idx_val)
            continue
            
        # 2. 3D spatial buffer exclusion
        if r['has_3d_coord'] == 1 and len(test_cas) > 0:
            ca_p = np.array([r['coord_x'], r['coord_y'], r['coord_z']])
            spatial_dist = np.min(np.linalg.norm(test_cas - ca_p, axis=1))
            if spatial_dist < spatial_buf:
                drop_tr_idx.append(idx_val)
                continue
                
    final_tr_idx = cand_tr.index.difference(drop_tr_idx).values
    return final_tr_idx, te_idx

target_configs = [
    ('UL30', 'Aciclovir', ['Exo', 'Palm', 'Fingers', 'Thumb_Cterm', 'NTD']),
    ('UL30', 'Foscarnet', ['Exo', 'Palm', 'Fingers', 'Thumb_Cterm', 'NTD']),
    ('UL23', 'Aciclovir', ['ATP_Core', 'Nucleoside_Pocket', 'Lid_Dimer', 'C-term', 'N-term'])
]

results = {}

for gene, drug, domain_list in target_configs:
    print(f"\nEvaluating {gene}-{drug} across {domain_list}...", flush=True)
    sub = df[(df['gene'] == gene) & (df['drug'] == drug)].copy()
    y_all = sub['target_binary'].values
    
    oof_pos_rf = np.full(len(sub), np.nan)
    oof_lr = np.full(len(sub), np.nan)
    oof_hgb = np.full(len(sub), np.nan)
    
    # Supervised feature set: PLM zero-shots + Where + Severity + Conservation + Radicalness
    sup_features = ['neg_esm2_650m', 'neg_protbert', 'neg_esm2_35m', 'weighted_conservation', 'blosum_radicalness', 'where_prereg', 'sev_prereg']
    for c in sup_features:
        if sub[c].isna().any():
            sub[c] = sub[c].fillna(sub[c].median())
            
    for d in domain_list:
        tr_idx, te_idx = create_dual_buffer_split(sub, d, seq_buf=5, spatial_buf=12.0)
        
        target_te_indices = sub.index.get_indexer(te_idx)
        
        X_pos_tr = sub.loc[tr_idx, ['pos']].values
        y_tr = sub.loc[tr_idx, 'target_binary'].values
        X_pos_te = sub.loc[te_idx, ['pos']].values
        
        n_tr_pos = int(y_tr.sum())
        n_tr_neg = int(len(y_tr) - n_tr_pos)
        print(f"  Domain {d}: Test N={len(te_idx)} (Pos={sub.loc[te_idx, 'target_binary'].sum():.0f}) | Train N={len(tr_idx)} (Pos={n_tr_pos}, Neg={n_tr_neg})", flush=True)
        
        if len(np.unique(y_tr)) > 1:
            rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5)
            rf.fit(X_pos_tr, y_tr)
            oof_pos_rf[target_te_indices] = rf.predict_proba(X_pos_te)[:, 1]
            
            # Supervised Logistic Regression L2
            X_tr_sup = sub.loc[tr_idx, sup_features].values
            X_te_sup = sub.loc[te_idx, sup_features].values
            
            mu_tr = np.mean(X_tr_sup, axis=0)
            std_tr = np.std(X_tr_sup, axis=0)
            std_tr[std_tr == 0] = 1.0
            
            X_tr_scaled = (X_tr_sup - mu_tr) / std_tr
            X_te_scaled = (X_te_sup - mu_tr) / std_tr
            
            lr = LogisticRegression(C=1.0, max_iter=1000, random_state=42)
            lr.fit(X_tr_scaled, y_tr)
            oof_lr[target_te_indices] = lr.predict_proba(X_te_scaled)[:, 1]
            
            hgb = HistGradientBoostingClassifier(max_iter=100, random_state=42, min_samples_leaf=5)
            hgb.fit(X_tr_sup, y_tr)
            oof_hgb[target_te_indices] = hgb.predict_proba(X_te_sup)[:, 1]
        else:
            base_rate = float(np.mean(y_tr)) if len(y_tr) > 0 else 0.5
            oof_pos_rf[target_te_indices] = base_rate
            oof_lr[target_te_indices] = base_rate
            oof_hgb[target_te_indices] = base_rate
            
    # Compute Pooled Out-Of-Fold AUROCs
    auc_pos_rf = roc_auc_score(y_all, oof_pos_rf)
    auc_lr = roc_auc_score(y_all, oof_lr)
    auc_hgb = roc_auc_score(y_all, oof_hgb)
    auc_zero_esm2 = roc_auc_score(y_all, sub['neg_esm2_650m'])
    auc_zero_protbert = roc_auc_score(y_all, sub['neg_protbert'])
    auc_zero_blosum = roc_auc_score(y_all, sub['blosum_radicalness'])
    
    results[f"{gene}-{drug}"] = {
        'n_total': int(len(sub)),
        'n_pos': int(y_all.sum()),
        'auc_pos_rf': float(auc_pos_rf),
        'auc_zero_esm2_650m': float(auc_zero_esm2),
        'auc_zero_protbert': float(auc_zero_protbert),
        'auc_zero_blosum': float(auc_zero_blosum),
        'auc_sup_logistic_regression': float(auc_lr),
        'auc_sup_hist_gradient_boosting': float(auc_hgb)
    }

print("\n" + "="*80, flush=True)
print("SUPERVISED VS ZERO-SHOT MODELS UNDER DUAL-BUFFERED LEAVE-ONE-DOMAIN-OUT", flush=True)
print("="*80, flush=True)
for k, v in results.items():
    print(f"\n--- Target: {k} (N={v['n_total']}, Pos={v['n_pos']}) ---", flush=True)
    print(f"Position-Only RF (Validity Check): {v['auc_pos_rf']:.4f}", flush=True)
    print(f"Zero-Shot ESM-2 650M:             {v['auc_zero_esm2_650m']:.4f}", flush=True)
    print(f"Zero-Shot ProtBERT:               {v['auc_zero_protbert']:.4f}", flush=True)
    print(f"Zero-Shot BLOSUM62 Radicalness:   {v['auc_zero_blosum']:.4f}", flush=True)
    print(f"Supervised Logistic Regression L2: {v['auc_sup_logistic_regression']:.4f}", flush=True)
    print(f"Supervised HistGradientBoosting:   {v['auc_sup_hist_gradient_boosting']:.4f}", flush=True)

with open('logs/supervised_dual_buffer_benchmark.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved results to logs/supervised_dual_buffer_benchmark.json", flush=True)
