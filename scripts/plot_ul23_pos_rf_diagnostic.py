"""
plot_ul23_pos_rf_diagnostic.py
Inspects position-only RF out-of-domain predictions vs position for each UL23 domain.
Saves plot to figures/ul23_pos_rf_diagnostic.png and prints numerical values.
"""
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

# 1. Load data
df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
sub = df[(df['gene'] == 'UL23') & (df['drug'] == 'Aciclovir')].copy().reset_index(drop=True)

evaluable_domains = ['ATP_Core', 'Nucleoside_Pocket', 'Lid_Dimer', 'C-term']

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

domain_diagnostics = {}

for idx, dom in enumerate(evaluable_domains):
    te_m = (sub['domain_5'] == dom)
    tr_m = (~te_m)
    
    # Dual buffer exclusion: seq_buf=5
    test_pos = sub.loc[te_m, 'pos'].values
    tr_drop = np.zeros(len(sub), dtype=bool)
    for p in test_pos:
        tr_drop |= (np.abs(sub['pos'].values - p) <= 5)
    tr_m = tr_m & (~tr_drop)
    
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

print("UL23 Position-Only RF Out-of-Domain Diagnostic:")
print("="*85)
for dom, d in domain_diagnostics.items():
    print(f"Domain: {dom:<18} | Range: {d['pos_range']} | N={d['n']} (Pos={d['pos']}) | AUROC = {d['auroc']:.4f}")
    print(f"  Mean Pred for Resistant: {d['pred_for_res']:.4f} | Mean Pred for Sensitive: {d['pred_for_sens']:.4f} | Diff: {d['pred_for_res'] - d['pred_for_sens']:+.4f}")

