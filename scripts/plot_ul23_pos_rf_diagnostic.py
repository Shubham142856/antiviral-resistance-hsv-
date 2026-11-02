"""
plot_ul23_pos_rf_diagnostic.py
Inspects position-only RF out-of-domain predictions vs position for each UL23 domain
under the exact dual-buffered protocol (seq_buf=5, spatial_buf=12.0 A).
Reconciles exactly with Table 7:
  - ATP_Core: 0.4013
  - Nucleoside_Pocket: 0.3815
  - Lid_Dimer: 0.4660
  - C-term: 0.5000
Saves plot to figures/ul23_pos_rf_diagnostic.png.
"""
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

# 1. Load data
df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
spatial_df = pd.read_csv('data/processed/primary_3d_spatial_folds.tsv', sep='\t')
df = pd.merge(df, spatial_df[['gene', 'pos', 'wt', 'mut', 'drug', 'coord_x', 'coord_y', 'coord_z', 'has_3d_coord']], on=['gene', 'pos', 'wt', 'mut', 'drug'], how='left')

sub = df[(df['gene'] == 'UL23') & (df['drug'] == 'Aciclovir')].copy().reset_index(drop=True)

evaluable_domains = ['ATP_Core', 'Nucleoside_Pocket', 'Lid_Dimer', 'C-term']

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

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

domain_diagnostics = {}

for idx, dom in enumerate(evaluable_domains):
    tr_m, te_m = get_dual_buffer_mask(sub, dom, seq_buf=5, spatial_buf=12.0)
    
    y_tr = sub.loc[tr_m, 'target_binary'].values
    y_te = sub.loc[te_m, 'target_binary'].values
    pos_te = sub.loc[te_m, 'pos'].values
    
    rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=1)
    rf.fit(sub.loc[tr_m, ['pos']], y_tr)
    p_pred = rf.predict_proba(sub.loc[te_m, ['pos']])[:, 1]
    auc = roc_auc_score(y_te, p_pred)
    
    # Sort by position for plotting
    s_idx = np.argsort(pos_te)
    pos_sorted = pos_te[s_idx]
    pred_sorted = p_pred[s_idx]
    y_sorted = y_te[s_idx]
    
    ax = axes[idx]
    ax.scatter(pos_sorted[y_sorted == 0], pred_sorted[y_sorted == 0], color='blue', alpha=0.7, label='Sensitive (0)')
    ax.scatter(pos_sorted[y_sorted == 1], pred_sorted[y_sorted == 1], color='red', alpha=0.7, label='Resistant (1)')
    ax.plot(pos_sorted, pred_sorted, color='black', linestyle='--', alpha=0.5)
    ax.axhline(0.5, color='gray', linestyle=':')
    ax.set_title(f"{dom} (N={len(y_te)}, Pos={int(y_te.sum())}, AUROC={auc:.4f})")
    ax.set_xlabel("Amino Acid Position")
    ax.set_ylabel("Pos-Only RF Predicted Probability")
    ax.legend(loc='best')
    ax.set_ylim(-0.05, 1.05)
    
    domain_diagnostics[dom] = {
        'n': len(y_te),
        'pos': int(y_te.sum()),
        'auroc': float(auc),
        'min_pred': float(np.min(p_pred)),
        'max_pred': float(np.max(p_pred)),
        'mean_pred': float(np.mean(p_pred)),
        'pred_for_res': float(np.mean(p_pred[y_te == 1])),
        'pred_for_sens': float(np.mean(p_pred[y_te == 0])),
        'pos_range': [int(np.min(pos_te)), int(np.max(pos_te))]
    }

plt.tight_layout()
pathlib.Path('figures').mkdir(exist_ok=True)
plt.savefig('figures/ul23_pos_rf_diagnostic.png', dpi=300)
plt.close()

print("UL23 Position-Only RF Out-of-Domain Diagnostic (Exact Dual-Buffer):")
print("="*85)
for dom, d in domain_diagnostics.items():
    print(f"Domain: {dom:<18} | Range: {d['pos_range']} | N={d['n']} (Pos={d['pos']}) | AUROC = {d['auroc']:.4f}")
    print(f"  Mean Pred for Resistant: {d['pred_for_res']:.4f} | Mean Pred for Sensitive: {d['pred_for_sens']:.4f} | Diff: {d['pred_for_res'] - d['pred_for_sens']:+.4f}")
