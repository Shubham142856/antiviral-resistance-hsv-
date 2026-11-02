"""
run_beta5_influential_positions_loo.py
Step 0b Item 1:
Investigate bimodal distribution of beta5 in UL30-ACV spec-exact model:
Model: logit = b0 + b1*W + b2*S + b3*S^2 + b4*(W*S) + b5*(W*S^2)
with penalty='l2', C=1.0 (exact match to logs/where_what_spec_exact_results.json).

1. Run Leave-One-Position-Out (LOO) for all 290 unique positions on UL30-ACV.
2. Identify top 15 most influential positions on beta5.
3. Quantify how presence/absence of these positions drives the two modes (~ -0.03 vs ~ 0.6-0.8).
4. Label: EXPLORATORY. Model and pre-registered falsifier remain unchanged.
"""
import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

# 1. Load data
df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
sub = df[(df['gene'] == 'UL30') & (df['drug'] == 'Aciclovir')].copy().reset_index(drop=True)

# 2. Spec-exact features
z_mu = - (sub['mu_plm_19'] - sub['mu_plm_19'].mean()) / sub['mu_plm_19'].std()
z_dist = (sub['dist_to_site_angstrom'] - sub['dist_to_site_angstrom'].mean()) / sub['dist_to_site_angstrom'].std()
sub['where_specexact'] = (z_mu - z_dist) / np.std(z_mu - z_dist)
sub['sev_specexact'] = - (sub['severity_plm'] - sub['severity_plm'].mean()) / sub['severity_plm'].std()

W = sub['where_specexact'].values
S = sub['sev_specexact'].values
y_full = sub['target_binary'].values.astype(int)

W_std = (W - np.mean(W)) / np.std(W)
S_std = (S - np.mean(S)) / np.std(S)

X_full = np.column_stack([
    W_std,
    S_std,
    S_std**2,
    W_std * S_std,
    W_std * (S_std**2)
])

# Fit baseline spec-exact model
lr_base = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42)
lr_base.fit(X_full, y_full)
base_beta5 = float(lr_base.coef_[0, 4])
base_coefs = [float(x) for x in lr_base.coef_[0]]

print(f"Baseline UL30-ACV beta5 (full cohort N={len(sub)}, Pos=62, UniquePos={len(sub['pos'].unique())}): {base_beta5:.4f}")

# Leave-One-Position-Out
unique_pos = np.sort(sub['pos'].unique())
loo_results = []

for p in unique_pos:
    mask = (sub['pos'] != p).values
    X_loo = X_full[mask]
    y_loo = y_full[mask]
    
    if len(np.unique(y_loo)) < 2:
        continue
        
    try:
        lr_loo = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42)
        lr_loo.fit(X_loo, y_loo)
        b5 = float(lr_loo.coef_[0, 4])
        shift = b5 - base_beta5
        
        p_sub = sub[sub['pos'] == p]
        dom = p_sub['domain_5'].iloc[0] if 'domain_5' in p_sub.columns else 'Unknown'
        wt = p_sub['wt'].iloc[0]
        n_var = len(p_sub)
        n_res = int(p_sub['target_binary'].sum())
        
        loo_results.append({
            'pos': int(p),
            'wt': str(wt),
            'domain': str(dom),
            'n_variants': int(n_var),
            'n_resistant': int(n_res),
            'beta5_loo': b5,
            'delta_beta5': shift,
            'abs_delta_beta5': abs(shift)
        })
    except Exception as e:
        continue

# Sort by absolute influence on beta5
loo_sorted = sorted(loo_results, key=lambda x: x['abs_delta_beta5'], reverse=True)
top15 = loo_sorted[:15]

print("\n" + "="*85)
print("TOP 15 INFLUENTIAL POSITIONS ON BETA5 (LEAVE-ONE-POSITION-OUT DIAGNOSTIC)")
print("Label: EXPLORATORY (Model and pre-registered falsifier unchanged)")
print("="*85)
print(f"{'Pos':<6} | {'WT':<4} | {'Domain':<12} | {'Variants (Res/Total)':<22} | {'Beta5 (LOO)':<12} | {'Delta Beta5':<12}")
print("-" * 85)
for r in top15:
    print(f"{r['pos']:<6} | {r['wt']:<4} | {r['domain']:<12} | {r['n_resistant']}/{r['n_variants']:<20} | {r['beta5_loo']:+.4f}     | {r['delta_beta5']:+.4f}")

# Bootstrap inclusion tracking
rng = np.random.RandomState(42)
pos_map = {p: np.where(sub['pos'].values == p)[0] for p in unique_pos}
n_boot = 2000
boot_b5 = []
top_pos_set = set([r['pos'] for r in top15[:5]])
top5_counts = []
resample_pos_inclusions = []

for _ in range(n_boot):
    s_p = rng.choice(unique_pos, size=len(unique_pos), replace=True)
    idx = np.concatenate([pos_map[p] for p in s_p])
    y_b = y_full[idx]
    if len(np.unique(y_b)) > 1:
        try:
            lr_b = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=500, random_state=42)
            lr_b.fit(X_full[idx], y_b)
            b = float(lr_b.coef_[0, 4])
            boot_b5.append(b)
            n_top_in = sum([np.sum(s_p == tp) for tp in top_pos_set])
            top5_counts.append(n_top_in)
            resample_pos_inclusions.append(s_p)
        except:
            pass

boot_b5 = np.array(boot_b5)
top5_counts = np.array(top5_counts)

# Bimodality threshold (inspection shows hump 1 around -0.03, hump 2 around +0.55-0.85)
mode1_mask = (boot_b5 < 0.3)
mode2_mask = (boot_b5 >= 0.3)

print("\n" + "="*85)
print("BIMODAL BOOTSTRAP DECOMPOSITION")
print("="*85)
print(f"Total valid resamples: {len(boot_b5)}")
print(f"Mode 1 (Beta5 < 0.3, centered ~ -0.03): {np.sum(mode1_mask)} ({np.mean(mode1_mask)*100:.1f}%) | Mean Beta5 = {np.mean(boot_b5[mode1_mask]):+.4f}")
print(f"Mode 2 (Beta5 >= 0.3, centered ~ +0.6 to +0.8): {np.sum(mode2_mask)} ({np.mean(mode2_mask)*100:.1f}%) | Mean Beta5 = {np.mean(boot_b5[mode2_mask]):+.4f}")
print(f"Mean inclusion of top-5 influential positions in Mode 1: {np.mean(top5_counts[mode1_mask]):.2f} copies")
print(f"Mean inclusion of top-5 influential positions in Mode 2: {np.mean(top5_counts[mode2_mask]):.2f} copies")
corr = np.corrcoef(boot_b5, top5_counts)[0, 1]
print(f"Correlation between top-5 influential positions count and Beta5: r = {corr:+.4f}")

# Check correlation of each of the top 15 positions individually with beta5
per_pos_corr = []
for r in top15:
    p = r['pos']
    p_counts = np.array([np.sum(s_p == p) for s_p in resample_pos_inclusions])
    r_corr = float(np.corrcoef(boot_b5, p_counts)[0, 1])
    per_pos_corr.append({
        'pos': p,
        'wt': r['wt'],
        'domain': r['domain'],
        'correlation_with_beta5': r_corr
    })
    print(f"  Pos {p} ({r['wt']}, {r['domain']}): correlation with beta5 in resamples = {r_corr:+.4f}")

out = {
    'label': 'EXPLORATORY',
    'status': 'Model and falsifier unchanged. Beta5 unidentifiable at current sample size.',
    'baseline_coefs': {
        'where_b1': base_coefs[0],
        'severity_b2': base_coefs[1],
        'severity_sq_b3': base_coefs[2],
        'where_x_severity_b4': base_coefs[3],
        'where_x_severity_sq_b5': base_coefs[4]
    },
    'top15_influential_positions': top15,
    'per_pos_resample_correlation': per_pos_corr,
    'bimodal_decomposition': {
        'total_resamples': int(len(boot_b5)),
        'mode1_count': int(np.sum(mode1_mask)),
        'mode1_pct': float(np.mean(mode1_mask)*100),
        'mode1_mean_beta5': float(np.mean(boot_b5[mode1_mask])),
        'mode2_count': int(np.sum(mode2_mask)),
        'mode2_pct': float(np.mean(mode2_mask)*100),
        'mode2_mean_beta5': float(np.mean(boot_b5[mode2_mask])),
        'correlation_top5_count_with_beta5': float(corr)
    }
}

with open('logs/beta5_top15_loo_diagnostic.json', 'w') as f:
    json.dump(out, f, indent=2)

print("\nSaved diagnostics to logs/beta5_top15_loo_diagnostic.json")
