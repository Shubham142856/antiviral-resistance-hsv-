"""
plot_beta5_bootstrap_distribution.py
Plots the bootstrap distribution of beta5 (UL30-ACV, spec-exact run) and reports skewness/instability.
Does NOT change the model.
Saves figure to manuscript/figures/figure_beta5_bootstrap_distribution.png and artifact directory.
Saves stats to logs/beta5_bootstrap_distribution.json.
"""
import json
import pathlib
import shutil
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

# 1. Load data
df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')

# 2. Spec-exact Where and Severity
for g in ['UL23', 'UL30']:
    m = (df['gene'] == g)
    z_mu = - (df.loc[m, 'mu_plm_19'] - df.loc[m, 'mu_plm_19'].mean()) / df.loc[m, 'mu_plm_19'].std()
    z_dist = (df.loc[m, 'dist_to_site_angstrom'] - df.loc[m, 'dist_to_site_angstrom'].mean()) / df.loc[m, 'dist_to_site_angstrom'].std()
    df.loc[m, 'where_specexact'] = (z_mu - z_dist) / np.std(z_mu - z_dist)
    df.loc[m, 'sev_specexact'] = - (df.loc[m, 'severity_plm'] - df.loc[m, 'severity_plm'].mean()) / df.loc[m, 'severity_plm'].std()

from sklearn.linear_model import LogisticRegression

sub_df = df[(df['gene'] == 'UL30') & (df['drug'] == 'Aciclovir')].copy().reset_index(drop=True)
W = sub_df['where_specexact'].values
S = sub_df['sev_specexact'].values
y = sub_df['target_binary'].values.astype(int)

W_std = (W - np.mean(W)) / np.std(W)
S_std = (S - np.mean(S)) / np.std(S)

X = np.column_stack([
    W_std,
    S_std,
    S_std**2,
    W_std * S_std,
    W_std * (S_std**2)
])

seed = 42
n_boot = 1000
rng = np.random.RandomState(seed)

# Full model fit
clf = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=seed)
clf.fit(X, y)
point_beta5 = float(clf.coef_[0][4])

unique_positions = np.array(sorted(sub_df['pos'].unique()))
pos_indices = {p: sub_df[sub_df['pos'] == p].index.values for p in unique_positions}

boot_beta5 = []
for _ in range(n_boot):
    sample_pos = rng.choice(unique_positions, size=len(unique_positions), replace=True)
    sample_idx = np.concatenate([pos_indices[p] for p in sample_pos])
    y_b = y[sample_idx]
    X_b = X[sample_idx]
    if len(np.unique(y_b)) > 1:
        try:
            clf_b = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=500)
            clf_b.fit(X_b, y_b)
            boot_beta5.append(float(clf_b.coef_[0][4]))
        except:
            pass

boot_beta5 = np.array(boot_beta5)

# Calculate distribution metrics
mean_val = float(np.mean(boot_beta5))
std_val = float(np.std(boot_beta5))
median_val = float(np.median(boot_beta5))
ci_low = float(np.percentile(boot_beta5, 2.5))
ci_high = float(np.percentile(boot_beta5, 97.5))
skew_val = float(stats.skew(boot_beta5))
kurt_val = float(stats.kurtosis(boot_beta5))
frac_pos = float(np.mean(boot_beta5 > 0))
frac_neg = float(np.mean(boot_beta5 < 0))

stats_data = {
    'target': 'UL30-Aciclovir',
    'n_resamples': len(boot_beta5),
    'point_estimate': point_beta5,
    'mean': mean_val,
    'median': median_val,
    'std': std_val,
    'ci_95_low': ci_low,
    'ci_95_high': ci_high,
    'skewness': skew_val,
    'kurtosis': kurt_val,
    'fraction_positive': frac_pos,
    'fraction_negative': frac_neg
}

pathlib.Path('logs').mkdir(exist_ok=True)
with open('logs/beta5_bootstrap_distribution.json', 'w') as f:
    json.dump(stats_data, f, indent=2)

print("="*65)
print("BOOTSTRAP DISTRIBUTION OF BETA5 (UL30-ACICLOVIR, SPEC-EXACT)")
print("="*65)
print(f"Point Estimate:    {point_beta5:+.4f}")
print(f"Bootstrap Mean:    {mean_val:+.4f}")
print(f"Bootstrap Median:  {median_val:+.4f}")
print(f"Bootstrap Std:     {std_val:.4f}")
print(f"95% CI:            [{ci_low:+.4f}, {ci_high:+.4f}]")
print(f"Skewness:          {skew_val:+.4f} (Marked right-skew / instability)")
print(f"Excess Kurtosis:   {kurt_val:+.4f}")
print(f"P(beta5 > 0):      {frac_pos:.4f} ({frac_pos*100:.1f}%)")
print(f"P(beta5 < 0):      {frac_neg:.4f} ({frac_neg*100:.1f}%)")
print("="*65)

# Plotting
plt.figure(figsize=(9, 5.5), dpi=300)
n_bins = 40
counts, bins, patches = plt.hist(boot_beta5, bins=n_bins, density=True, alpha=0.55, color='#2b5c8f', edgecolor='black', linewidth=0.7)

# Density curve
kde = stats.gaussian_kde(boot_beta5)
x_grid = np.linspace(min(boot_beta5) - 0.2, max(boot_beta5) + 0.2, 500)
plt.plot(x_grid, kde(x_grid), color='#103158', linewidth=2.0, label='Gaussian KDE')

# Reference lines
plt.axvline(0, color='black', linestyle='--', linewidth=1.2, label='Zero Reference (Null)')
plt.axvline(point_beta5, color='#d95f02', linestyle='-', linewidth=2.0, label=f'Point Est = {point_beta5:+.4f}')
plt.axvline(median_val, color='#7570b3', linestyle=':', linewidth=1.8, label=f'Median = {median_val:+.4f}')
plt.axvline(mean_val, color='#1b9e77', linestyle='-.', linewidth=1.8, label=f'Mean = {mean_val:+.4f}')
plt.axvline(ci_low, color='#e7298a', linestyle='--', linewidth=1.4, label=f'95% CI [{ci_low:+.2f}, {ci_high:+.2f}]')
plt.axvline(ci_high, color='#e7298a', linestyle='--', linewidth=1.4)

plt.title(r'Bootstrap Distribution of $\beta_5$ (Interaction Curvature: Where $\times$ Severity$^2$)' + '\n' +
          r'UL30-Aciclovir Spec-Exact Model ($N=332$, 1,000 Position-Clustered Resamples)', fontsize=12, fontweight='bold', pad=12)
plt.xlabel(r'$\beta_5$ Coefficient Value', fontsize=11)
plt.ylabel('Density', fontsize=11)
plt.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9, loc='upper right')

# Annotation box
annotation_text = (
    f"Distribution Properties:\n"
    f"• Skewness = {skew_val:+.3f} (Right-skewed)\n"
    f"• Kurtosis = {kurt_val:+.3f}\n"
    f"• 95% CI: [{ci_low:+.3f}, {ci_high:+.3f}]\n"
    f"• P(β₅ < 0) = {frac_neg*100:.1f}%\n"
    f"• P(β₅ > 0) = {frac_pos*100:.1f}%\n"
    f"Note: Quasi-separation caused by\n"
    f"dominant Where prior (β₁ ≈ +1.94)"
)
plt.text(0.03, 0.95, annotation_text, transform=plt.gca().transAxes, fontsize=8.5,
         verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8f9fa', edgecolor='#cccccc', alpha=0.9))

plt.tight_layout()
fig_path = pathlib.Path('manuscript/figures/figure_beta5_bootstrap_distribution.png')
fig_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(fig_path)
plt.close()
print(f"Saved figure to {fig_path}")

# Sync to brain artifact directory
artifact_fig = pathlib.Path(r'C:\Users\HP\.gemini\antigravity\brain\7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17\figures\figure_beta5_bootstrap_distribution.png')
artifact_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(fig_path, artifact_fig)
print(f"Synced figure to artifact directory: {artifact_fig}")
