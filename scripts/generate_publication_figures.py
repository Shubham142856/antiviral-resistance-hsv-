"""
Generate Publication Figures 1-4 for HSV Antiviral Resistance Benchmark
========================================================================
Visualizes Tables 7-11 and modern confirmatory / exploratory findings at 300 DPI:
  - Figure 1: Spatial Shortcut Collapse & Sequence-Neighbour Gap Metric (Tables 11A, 11B)
  - Figure 2: The Evidence Ladder & PSSM Ascertainment Robustness (Table 8)
  - Figure 3: Confirmatory Exonuclease Double Dissociation & Leave-Site-Out Gradient (Tables 9E, 9F)
  - Figure 4: Feature-Family Degradation Matrix & Out-of-Block Generalization (Table 11B, Paired Contrasts)

Reads strictly verified, committed JSON logs:
  - logs/split_separation_diagnostic.json
  - logs/split_degradation_benchmark_results.json
  - logs/evidence_ladder_table8_results.json
  - logs/step2_hcmv_exonuclease_confirmatory_results.json
  - logs/ul54_leave_site_out_robustness.json
  - logs/paired_family_contrast_results.json

Outputs high-DPI publication figures to:
  1. Artifacts directory: C:/Users/HP/.gemini/antigravity/brain/7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17/figures/
  2. Workspace directory: manuscript/figures/
"""

import os
import json
import shutil
import numpy as np
import matplotlib.pyplot as plt

# Publication aesthetic styling (Nature / ICML standard)
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['grid.color'] = '#E0E0E0'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.linewidth'] = 0.5
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.framealpha'] = 0.9
plt.rcParams['legend.edgecolor'] = '#CCCCCC'

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LOGS = os.path.join(ROOT, "logs")

ARTIFACT_DIR = r"C:\Users\HP\.gemini\antigravity\brain\7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17\figures"
WORKSPACE_DIR = os.path.join(ROOT, "manuscript", "figures")

os.makedirs(ARTIFACT_DIR, exist_ok=True)
os.makedirs(WORKSPACE_DIR, exist_ok=True)

def save_fig(fig, base_name, legacy_names=None):
    for d in [ARTIFACT_DIR, WORKSPACE_DIR]:
        path = os.path.join(d, base_name)
        fig.savefig(path, bbox_inches='tight', dpi=300)
        print(f"Saved: {path}")
        if legacy_names:
            for leg in legacy_names:
                leg_path = os.path.join(d, leg)
                shutil.copyfile(path, leg_path)
                print(f"Saved legacy copy: {leg_path}")
    plt.close(fig)

# =============================================================================
# FIGURE 1: Spatial Shortcut Collapse & Sequence-Neighbour Gap Metric
# =============================================================================
print("Generating Figure 1: Spatial Shortcut Collapse & Sequence Gap Metric...")
with open(os.path.join(LOGS, "split_separation_diagnostic.json"), "r") as f:
    diag_data = json.load(f)
with open(os.path.join(LOGS, "split_degradation_benchmark_results.json"), "r") as f:
    degrad_data = json.load(f)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0), dpi=300)

# Panel A: Sequence-distance diagnostic (Table 11A)
targets_order = ['UL23-Aciclovir', 'UL30-Aciclovir', 'UL30-Foscarnet']
target_labels = ['UL23-ACV\n(Thymidine Kinase)', 'UL30-ACV\n(DNA Polymerase)', 'UL30-FOS\n(DNA Polymerase)']

pg_medians = [diag_data['targets'][t]['schemes']['position-grouped']['median_nearest_train_distance'] for t in targets_order]
cb_medians = [diag_data['targets'][t]['schemes']['contiguous-block']['median_nearest_train_distance'] for t in targets_order]
pg_frac2 = [diag_data['targets'][t]['schemes']['position-grouped']['frac_within_2_residues'] * 100 for t in targets_order]
cb_frac2 = [diag_data['targets'][t]['schemes']['contiguous-block']['frac_within_2_residues'] * 100 for t in targets_order]

x = np.arange(len(targets_order))
w = 0.35

rects1 = ax1.bar(x - w/2, pg_medians, w, label='Position-Grouped CV (Standard)',
                 color='#E65100', edgecolor='#BF360C', alpha=0.9)
rects2 = ax1.bar(x + w/2, cb_medians, w, label='Contiguous-Block Holdout (Isolating)',
                 color='#1976D2', edgecolor='#0D47A1', alpha=0.9)

ax1.set_ylabel('Median Nearest Training Distance (Residues)', fontsize=10.5, fontweight='bold')
ax1.set_title('A. Sequence-Neighbour Separation Diagnostic (Table 11A)', fontsize=11.5, fontweight='bold', pad=10)
ax1.set_xticks(x)
ax1.set_xticklabels(target_labels, fontsize=9.5)
ax1.set_yscale('log')
ax1.set_ylim(0.5, 200)
ax1.grid(True, axis='y', alpha=0.6)
ax1.legend(loc='upper left', fontsize=9)

# Value annotations with % within 2 aa
for i in range(len(targets_order)):
    ax1.text(x[i] - w/2, pg_medians[i] * 1.25, f"{pg_medians[i]:.1f} aa\n({pg_frac2[i]:.1f}% ≤2aa)",
             ha='center', va='bottom', fontsize=8, fontweight='bold', color='#BF360C')
    ax1.text(x[i] + w/2, cb_medians[i] * 1.15, f"{cb_medians[i]:.1f} aa\n({cb_frac2[i]:.1f}% ≤2aa)",
             ha='center', va='bottom', fontsize=8, fontweight='bold', color='#0D47A1')

ax1.annotate('Memorization Corridor:\n96.5% of UL23 test sites\nsit ≤2 aa from training',
             xy=(0 - w/2, 1.0), xytext=(-0.35, 12),
             arrowprops=dict(facecolor='#BF360C', arrowstyle='->', lw=1.2),
             fontsize=8.2, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#FFF3E0", ec="#E65100"))

# Panel B: Positional Identity AUROC Collapse (Table 11B)
# Compare Random Forest vs Logistic Regression on Positional Identity
pairs = ['UL30-ACV (RF)', 'UL30-FOS (RF)', 'UL23-ACV (RF)', 'UL30-ACV (LR)', 'UL30-FOS (LR)', 'UL23-ACV (LR)']
pg_pos_auc = [
    degrad_data['targets']['UL30-Aciclovir']['models']['random_forest']['Positional Identity']['position_grouped_auc'],
    degrad_data['targets']['UL30-Foscarnet']['models']['random_forest']['Positional Identity']['position_grouped_auc'],
    degrad_data['targets']['UL23-Aciclovir']['models']['random_forest']['Positional Identity']['position_grouped_auc'],
    degrad_data['targets']['UL30-Aciclovir']['models']['logistic_regression']['Positional Identity']['position_grouped_auc'],
    degrad_data['targets']['UL30-Foscarnet']['models']['logistic_regression']['Positional Identity']['position_grouped_auc'],
    degrad_data['targets']['UL23-Aciclovir']['models']['logistic_regression']['Positional Identity']['position_grouped_auc'],
]
cb_pos_auc = [
    degrad_data['targets']['UL30-Aciclovir']['models']['random_forest']['Positional Identity']['contiguous_block_auc'],
    degrad_data['targets']['UL30-Foscarnet']['models']['random_forest']['Positional Identity']['contiguous_block_auc'],
    degrad_data['targets']['UL23-Aciclovir']['models']['random_forest']['Positional Identity']['contiguous_block_auc'],
    degrad_data['targets']['UL30-Aciclovir']['models']['logistic_regression']['Positional Identity']['contiguous_block_auc'],
    degrad_data['targets']['UL30-Foscarnet']['models']['logistic_regression']['Positional Identity']['contiguous_block_auc'],
    degrad_data['targets']['UL23-Aciclovir']['models']['logistic_regression']['Positional Identity']['contiguous_block_auc'],
]

x2 = np.arange(len(pairs))
w2 = 0.38

rects3 = ax2.bar(x2 - w2/2, pg_pos_auc, w2, label='Position-Grouped CV (Spurious)',
                 color='#D32F2F', edgecolor='#B71C1C', alpha=0.88)
rects4 = ax2.bar(x2 + w2/2, cb_pos_auc, w2, label='Contiguous-Block Holdout (Isolated)',
                 color='#78909C', edgecolor='#37474F', alpha=0.88)

ax2.axhline(0.50, color='black', linestyle=':', linewidth=1.1, label='Chance Baseline (0.50)')
ax2.set_ylabel('Positional Identity AUROC', fontsize=10.5, fontweight='bold')
ax2.set_title('B. Positional Identity Shortcut Collapse (Table 11B)', fontsize=11.5, fontweight='bold', pad=10)
ax2.set_xticks(x2)
ax2.set_xticklabels(pairs, fontsize=8.5, rotation=15, ha='right')
ax2.set_ylim(0.15, 1.0)
ax2.grid(True, axis='y', alpha=0.6)
ax2.legend(loc='upper right', fontsize=8.5)

# Annotations for deltas
for i in range(len(pairs)):
    delta = pg_pos_auc[i] - cb_pos_auc[i]
    ax2.text(x2[i], max(pg_pos_auc[i], cb_pos_auc[i]) + 0.03, f"Δ={delta:+.2f}",
             ha='center', fontsize=7.8, fontweight='bold', color='#B71C1C')

ax2.annotate('Linear model inverts\nbelow chance (0.22-0.35)\ndue to extrapolation',
             xy=(3 + w2/2, cb_pos_auc[3]), xytext=(2.6, 0.32),
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.0),
             fontsize=8.0, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#ECEFF1", ec="#455A64"))

plt.tight_layout()
save_fig(fig, "figure1_spatial_shortcut_collapse.png", ["figure1_spatial_leakage_collapse.png"])

# =============================================================================
# FIGURE 2: The Evidence Ladder & PSSM Ascertainment Robustness
# =============================================================================
print("Generating Figure 2: Evidence Ladder & PSSM Ascertainment Robustness...")
with open(os.path.join(LOGS, "evidence_ladder_table8_results.json"), "r") as f:
    ladder_data = json.load(f)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0), dpi=300)

# Panel A: Ladder Trajectories across D1 -> D2 -> D3 -> S2
ladder_tiers = ['D1 Pooled', 'D2 Domain-Mean', 'D3 Within-Study (S0)', 'S2 IC50 Stratum']
ul30_pssm = [
    ladder_data['UL30-Aciclovir']['scores']['PSSM']['D1_pooled']['auc'],
    ladder_data['UL30-Aciclovir']['scores']['PSSM']['D2_domain_mean']['auc'],
    ladder_data['UL30-Aciclovir']['scores']['PSSM']['D3_within_study']['pair_weighted'],
    ladder_data['UL30-Aciclovir']['scores']['PSSM']['S2_numeric_fc']['auc']
]
ul30_esm = [
    ladder_data['UL30-Aciclovir']['scores']['ESM-2 650M']['D1_pooled']['auc'],
    ladder_data['UL30-Aciclovir']['scores']['ESM-2 650M']['D2_domain_mean']['auc'],
    ladder_data['UL30-Aciclovir']['scores']['ESM-2 650M']['D3_within_study']['pair_weighted'],
    ladder_data['UL30-Aciclovir']['scores']['ESM-2 650M']['S2_numeric_fc']['auc']
]
ul30_ent = [
    ladder_data['UL30-Aciclovir']['scores']['-entropy']['D1_pooled']['auc'],
    ladder_data['UL30-Aciclovir']['scores']['-entropy']['D2_domain_mean']['auc'],
    ladder_data['UL30-Aciclovir']['scores']['-entropy']['D3_within_study']['pair_weighted'],
    ladder_data['UL30-Aciclovir']['scores']['-entropy']['S2_numeric_fc']['auc']
]
ul30_blosum = [
    ladder_data['UL30-Aciclovir']['scores']['BLOSUM62']['D1_pooled']['auc'],
    ladder_data['UL30-Aciclovir']['scores']['BLOSUM62']['D2_domain_mean']['auc'],
    ladder_data['UL30-Aciclovir']['scores']['BLOSUM62']['D3_within_study']['pair_weighted'],
    ladder_data['UL30-Aciclovir']['scores']['BLOSUM62']['S2_numeric_fc']['auc']
]

ax1.plot(ladder_tiers, ul30_esm, marker='o', lw=2.2, color='#1565C0', label='ESM-2 650M (Zero-Shot)')
ax1.plot(ladder_tiers, ul30_pssm, marker='s', lw=2.2, color='#2E7D32', label='Multi-Herpesvirus PSSM')
ax1.plot(ladder_tiers, ul30_ent, marker='^', lw=1.8, color='#E65100', label='Positional Shannon Entropy')
ax1.plot(ladder_tiers, ul30_blosum, marker='x', lw=1.6, ls='--', color='#757575', label='BLOSUM62 (Negative Control)')
ax1.axhline(0.50, color='black', linestyle=':', linewidth=1.1, label='Chance (0.50)')

ax1.set_ylabel('Harmonized AUROC', fontsize=10.5, fontweight='bold')
ax1.set_title('A. Evidence Ladder Trajectory (UL30-Aciclovir, Table 8)', fontsize=11.5, fontweight='bold', pad=10)
ax1.set_ylim(0.40, 0.90)
ax1.grid(True, alpha=0.6)
ax1.legend(loc='lower left', fontsize=8.5)

ax1.annotate('Classical PSSM stabilizes at 0.7079\nstrictly matching ESM-2 (0.7079)',
             xy=(2, 0.7079), xytext=(1.2, 0.78),
             arrowprops=dict(facecolor='#2E7D32', arrowstyle='->', lw=1.2),
             fontsize=8.2, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#E8F5E9", ec="#2E7D32"))

# Panel B: PSSM vs ESM-2 650M Confidence Intervals Under D3 Within-Study S0
ladder_arms = ['UL30-ACV', 'UL30-FOS', 'UL23-ACV']
arm_keys = ['UL30-Aciclovir', 'UL30-Foscarnet', 'UL23-Aciclovir']

pssm_pts = [ladder_data[k]['scores']['PSSM']['D3_within_study']['pair_weighted'] for k in arm_keys]
pssm_cis = [ladder_data[k]['scores']['PSSM']['D3_within_study']['ci_95'] for k in arm_keys]
pssm_err_lo = [pt - ci[0] for pt, ci in zip(pssm_pts, pssm_cis)]
pssm_err_hi = [ci[1] - pt for pt, ci in zip(pssm_pts, pssm_cis)]

esm_pts = [ladder_data[k]['scores']['ESM-2 650M']['D3_within_study']['pair_weighted'] for k in arm_keys]
esm_cis = [ladder_data[k]['scores']['ESM-2 650M']['D3_within_study']['ci_95'] for k in arm_keys]
esm_err_lo = [pt - ci[0] for pt, ci in zip(esm_pts, esm_cis)]
esm_err_hi = [ci[1] - pt for pt, ci in zip(esm_pts, esm_cis)]

blosum_pts = [ladder_data[k]['scores']['BLOSUM62']['D3_within_study']['pair_weighted'] for k in arm_keys]
blosum_cis = [ladder_data[k]['scores']['BLOSUM62']['D3_within_study']['ci_95'] for k in arm_keys]
blosum_err_lo = [pt - ci[0] for pt, ci in zip(blosum_pts, blosum_cis)]
blosum_err_hi = [ci[1] - pt for pt, ci in zip(blosum_pts, blosum_cis)]

x_b = np.arange(len(ladder_arms))
w_b = 0.26

r1 = ax2.bar(x_b - w_b, pssm_pts, w_b, yerr=[pssm_err_lo, pssm_err_hi], capsize=4,
             label='PSSM Profile (3/3 Excludes 0.50)', color='#2E7D32', edgecolor='#1B5E20', alpha=0.9)
r2 = ax2.bar(x_b, esm_pts, w_b, yerr=[esm_err_lo, esm_err_hi], capsize=4,
             label='ESM-2 650M (Interval Widening)', color='#1565C0', edgecolor='#0D47A1', alpha=0.9)
r3 = ax2.bar(x_b + w_b, blosum_pts, w_b, yerr=[blosum_err_lo, blosum_err_hi], capsize=4,
             label='BLOSUM62 (Negative Control)', color='#9E9E9E', edgecolor='#616161', alpha=0.85)

ax2.axhline(0.50, color='red', linestyle='--', linewidth=1.1, label='Chance (0.50)')
ax2.set_ylabel('D3 Within-Study Pair-Weighted AUROC (95% CI)', fontsize=10.5, fontweight='bold')
ax2.set_title('B. Classical PSSM Robustness vs Deep Learning (D3 Tier)', fontsize=11.5, fontweight='bold', pad=10)
ax2.set_xticks(x_b)
ax2.set_xticklabels(ladder_arms, fontsize=9.5)
ax2.set_ylim(0.15, 1.05)
ax2.grid(True, axis='y', alpha=0.6)
ax2.legend(loc='lower left', fontsize=8.2)

# Annotation for 3/3 robustness
ax2.annotate('PSSM strictly excludes 0.50\nacross 3/3 targets\n[0.61-0.80, 0.62-0.93, 0.60-0.80]',
             xy=(0 - w_b, pssm_pts[0]), xytext=(-0.35, 0.88),
             arrowprops=dict(facecolor='#2E7D32', arrowstyle='->', lw=1.2),
             fontsize=8.0, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#E8F5E9", ec="#2E7D32"))
ax2.annotate('ESM-2 95% CI spans chance\non UL30-FOS [0.35, 0.90]',
             xy=(1, esm_pts[1]), xytext=(0.85, 0.24),
             arrowprops=dict(facecolor='#1565C0', arrowstyle='->', lw=1.2),
             fontsize=8.0, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#E3F2FD", ec="#1565C0"))

plt.tight_layout()
save_fig(fig, "figure2_evidence_ladder_pssm_robustness.png", ["figure2_biological_dissociation.png"])

# =============================================================================
# FIGURE 3: Confirmatory Exonuclease Double Dissociation & Leave-Site-Out Gradient
# =============================================================================
print("Generating Figure 3: Exonuclease Double Dissociation & Leave-Site-Out...")
with open(os.path.join(LOGS, "step2_hcmv_exonuclease_confirmatory_results.json"), "r") as f:
    exo_data = json.load(f)
with open(os.path.join(LOGS, "ul54_leave_site_out_robustness.json"), "r") as f:
    lso_data = json.load(f)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0), dpi=300)

# Panel A: Active-Site Proximity Crossover (Table 9E)
drugs = ['Ganciclovir', 'Cidofovir', 'Foscarnet']
drug_labels = ['UL54-Ganciclovir\n(Chain Terminator)', 'UL54-Cidofovir\n(Chain Terminator)', 'UL54-Foscarnet\n(Pyrophosphate Inh.)']

pol_pts = [exo_data['targets'][d]['s0_polymerase']['point'] for d in drugs]
pol_cis = [exo_data['targets'][d]['s0_polymerase']['ci_95'] for d in drugs]
pol_err_lo = [p - c[0] for p, c in zip(pol_pts, pol_cis)]
pol_err_hi = [c[1] - p for p, c in zip(pol_pts, pol_cis)]

exo_pts = [exo_data['targets'][d]['s0_exonuclease']['point'] for d in drugs]
exo_cis = [exo_data['targets'][d]['s0_exonuclease']['ci_95'] for d in drugs]
exo_err_lo = [p - c[0] for p, c in zip(exo_pts, exo_cis)]
exo_err_hi = [c[1] - p for p, c in zip(exo_pts, exo_cis)]

x_d = np.arange(len(drugs))
w_d = 0.35

rects_pol = ax1.bar(x_d - w_d/2, pol_pts, w_d, yerr=[pol_err_lo, pol_err_hi], capsize=4,
                    label='Polymerase Active-Site Distance', color='#D84315', edgecolor='#BF360C', alpha=0.9)
rects_exo = ax1.bar(x_d + w_d/2, exo_pts, w_d, yerr=[exo_err_lo, exo_err_hi], capsize=4,
                    label='Exonuclease Pocket Distance', color='#00897B', edgecolor='#004D40', alpha=0.9)

ax1.axhline(0.50, color='black', linestyle=':', linewidth=1.1, label='Chance (0.50)')
ax1.set_ylabel('S0 Pair-Weighted AUROC (95% CI)', fontsize=10.5, fontweight='bold')
ax1.set_title('A. Confirmatory Double Dissociation Crossover (Table 9E)', fontsize=11.5, fontweight='bold', pad=10)
ax1.set_xticks(x_d)
ax1.set_xticklabels(drug_labels, fontsize=9.2)
ax1.set_ylim(0.05, 1.05)
ax1.grid(True, axis='y', alpha=0.6)
ax1.legend(loc='lower left', fontsize=8.8)

# Paired deltas annotations
for i, d in enumerate(drugs):
    delta = exo_data['targets'][d]['delta_exo_minus_pol']['observed']
    p_v = exo_data['targets'][d]['delta_exo_minus_pol']['p_value']
    p_str = "p<0.0005" if p_v < 0.0005 else f"p={p_v:.3f}"
    ax1.text(x_d[i], max(pol_pts[i], exo_pts[i]) + 0.08, f"Δ = {delta:+.3f}\n({p_str})",
             ha='center', fontsize=8.0, fontweight='bold', color='#004D40' if delta > 0 else '#BF360C')

ax1.annotate('Exonuclease distance recovers\nchain terminators (+0.31 to +0.68)\nwhile inverting on FOS (-0.56)',
             xy=(1 + w_d/2, exo_pts[1]), xytext=(0.4, 0.22),
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.1),
             fontsize=8.0, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#E0F2F1", ec="#00897B"))

# Panel B: Leave-Site-Out Robustness (Table 9F)
# Compares S0 Full Cohort vs S0 Own Defining Residues Removed
lso_drugs = ['Ganciclovir', 'Cidofovir', 'Foscarnet']
lso_labels = ['UL54-GCV\n(Exo Site)', 'UL54-CDV\n(Exo Site)', 'UL54-FOS\n(Pol Site)']

full_s0 = []
drop_s0 = []
drop_cis = []

for d in lso_drugs:
    feat_key = lso_data['survival_verdicts'][d]['primary_feature']
    ci_key = feat_key + '_ci_95'
    full_s0.append(lso_data['drugs'][d]['cohorts']['full'][feat_key])
    drop_s0.append(lso_data['drugs'][d]['cohorts']['drop_both'][feat_key])
    drop_cis.append(lso_data['drugs'][d]['cohorts']['drop_both'][ci_key])

drop_err_lo = [p - c[0] for p, c in zip(drop_s0, drop_cis)]
drop_err_hi = [c[1] - p for p, c in zip(drop_s0, drop_cis)]

x_l = np.arange(len(lso_drugs))
w_l = 0.35

rects_full = ax2.bar(x_l - w_l/2, full_s0, w_l, label='Full Cohort S0',
                     color='#5C6BC0', edgecolor='#3949AB', alpha=0.88)
rects_drop = ax2.bar(x_l + w_l/2, drop_s0, w_l, yerr=[drop_err_lo, drop_err_hi], capsize=4,
                     label='Own Defining Residues Removed (Leave-Site-Out)',
                     color='#26A69A', edgecolor='#00695C', alpha=0.9)

ax2.axhline(0.50, color='red', linestyle='--', linewidth=1.1, label='Chance (0.50)')
ax2.axhline(0.55, color='gray', linestyle=':', linewidth=1.0, label='Survival Threshold (0.55)')
ax2.set_ylabel('S0 Pair-Weighted AUROC (95% CI)', fontsize=10.5, fontweight='bold')
ax2.set_title('B. Leave-Site-Out Spatial Gradient Robustness (Table 9F)', fontsize=11.5, fontweight='bold', pad=10)
ax2.set_xticks(x_l)
ax2.set_xticklabels(lso_labels, fontsize=9.2)
ax2.set_ylim(0.40, 1.05)
ax2.grid(True, axis='y', alpha=0.6)
ax2.legend(loc='lower left', fontsize=8.5)

for i in range(len(lso_drugs)):
    diff = drop_s0[i] - full_s0[i]
    ax2.text(x_l[i] + w_l/2, drop_s0[i] + 0.05, f"{drop_s0[i]:.4f}\n(Δ={diff:+.3f})",
             ha='center', fontsize=8.0, fontweight='bold', color='#00695C')

ax2.annotate('All 3 arms retain signal\n(Shift ≤ 0.045 AUROC;\nall CIs strictly exclude 0.50)',
             xy=(1 + w_l/2, drop_s0[1]), xytext=(0.6, 0.60),
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.1),
             fontsize=8.0, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#E0F2F1", ec="#00695C"))

plt.tight_layout()
save_fig(fig, "figure3_exonuclease_double_dissociation.png", ["figure3_model_benchmark_matrix.png"])

# =============================================================================
# FIGURE 4: Feature-Family Degradation Matrix & Out-of-Block Generalization
# =============================================================================
print("Generating Figure 4: Feature Degradation Matrix & Paired Contrasts...")
with open(os.path.join(LOGS, "paired_family_contrast_results.json"), "r") as f:
    pfc_data = json.load(f)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0), dpi=300)

# Panel A: Degradation Deltas (PG - Block) Across Families (Table 11B)
arms = ['UL30-ACV', 'UL30-FOS', 'UL23-ACV']
arm_keys = ['UL30-Aciclovir', 'UL30-Foscarnet', 'UL23-Aciclovir']

pos_deltas_rf = [degrad_data['targets'][k]['models']['random_forest']['Positional Identity']['degradation_delta_pg_minus_block'] for k in arm_keys]
cons_deltas_rf = [degrad_data['targets'][k]['models']['random_forest']['Evolutionary Conservation']['degradation_delta_pg_minus_block'] for k in arm_keys]
plm_deltas_rf = [degrad_data['targets'][k]['models']['random_forest']['Protein Language Models']['degradation_delta_pg_minus_block'] for k in arm_keys]

x_a = np.arange(len(arms))
w_a = 0.25

ax1.bar(x_a - w_a, pos_deltas_rf, w_a, label='Positional Identity (Severe Shortcut)',
        color='#D32F2F', edgecolor='#B71C1C', alpha=0.9)
ax1.bar(x_a, cons_deltas_rf, w_a, label='Evolutionary Conservation (Robust)',
        color='#1976D2', edgecolor='#0D47A1', alpha=0.9)
ax1.bar(x_a + w_a, plm_deltas_rf, w_a, label='Protein Language Models (Robust)',
        color='#388E3C', edgecolor='#1B5E20', alpha=0.9)

ax1.axhline(0.0, color='black', linestyle='-', linewidth=0.8)
ax1.set_ylabel('Degradation Delta: AUROC(PG) - AUROC(Block)', fontsize=10.5, fontweight='bold')
ax1.set_title('A. Supervised Feature-Family Degradation Matrix (RF, Table 11B)', fontsize=11.5, fontweight='bold', pad=10)
ax1.set_xticks(x_a)
ax1.set_xticklabels(arms, fontsize=9.5)
ax1.set_ylim(-0.05, 0.50)
ax1.grid(True, axis='y', alpha=0.6)
ax1.legend(loc='upper right', fontsize=8.5)

for i in range(len(arms)):
    ax1.text(x_a[i] - w_a, pos_deltas_rf[i] + 0.015, f"+{pos_deltas_rf[i]:.3f}", ha='center', fontsize=7.8, fontweight='bold', color='#B71C1C')
    ax1.text(x_a[i], cons_deltas_rf[i] + 0.015, f"+{cons_deltas_rf[i]:.3f}", ha='center', fontsize=7.8, fontweight='bold', color='#0D47A1')
    ax1.text(x_a[i] + w_a, plm_deltas_rf[i] + 0.015, f"+{plm_deltas_rf[i]:.3f}", ha='center', fontsize=7.8, fontweight='bold', color='#1B5E20')

ax1.annotate('Leakage resides in\nfeature representation:\nPositional: +0.27 to +0.42\nPLM & Conserv: ≤ +0.07',
             xy=(0 - w_a, pos_deltas_rf[0]), xytext=(-0.35, 0.32),
             arrowprops=dict(facecolor='#B71C1C', arrowstyle='->', lw=1.1),
             fontsize=8.0, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#FFEBEE", ec="#D32F2F"))

# Panel B: Out-of-Block Paired Contrasts (PLM minus Conservation)
contrast_pairs = [
    'UL30-ACV\n(Random Forest)',
    'UL30-ACV\n(Logistic Reg.)',
    'UL30-FOS\n(Random Forest)',
    'UL30-FOS\n(Logistic Reg.)',
    'UL23-ACV\n(Random Forest)',
    'UL23-ACV\n(Logistic Reg.)'
]

plm_minus_cons = [
    pfc_data['targets']['UL30-Aciclovir']['models']['random_forest']['paired_block_auroc_plm_minus_conservation']['mean'],
    pfc_data['targets']['UL30-Aciclovir']['models']['logistic_regression']['paired_block_auroc_plm_minus_conservation']['mean'],
    pfc_data['targets']['UL30-Foscarnet']['models']['random_forest']['paired_block_auroc_plm_minus_conservation']['mean'],
    pfc_data['targets']['UL30-Foscarnet']['models']['logistic_regression']['paired_block_auroc_plm_minus_conservation']['mean'],
    pfc_data['targets']['UL23-Aciclovir']['models']['random_forest']['paired_block_auroc_plm_minus_conservation']['mean'],
    pfc_data['targets']['UL23-Aciclovir']['models']['logistic_regression']['paired_block_auroc_plm_minus_conservation']['mean'],
]
plm_minus_cons_ci = [
    pfc_data['targets']['UL30-Aciclovir']['models']['random_forest']['paired_block_auroc_plm_minus_conservation']['ci_95'],
    pfc_data['targets']['UL30-Aciclovir']['models']['logistic_regression']['paired_block_auroc_plm_minus_conservation']['ci_95'],
    pfc_data['targets']['UL30-Foscarnet']['models']['random_forest']['paired_block_auroc_plm_minus_conservation']['ci_95'],
    pfc_data['targets']['UL30-Foscarnet']['models']['logistic_regression']['paired_block_auroc_plm_minus_conservation']['ci_95'],
    pfc_data['targets']['UL23-Aciclovir']['models']['random_forest']['paired_block_auroc_plm_minus_conservation']['ci_95'],
    pfc_data['targets']['UL23-Aciclovir']['models']['logistic_regression']['paired_block_auroc_plm_minus_conservation']['ci_95'],
]

err_lo_c = [m - c[0] for m, c in zip(plm_minus_cons, plm_minus_cons_ci)]
err_hi_c = [c[1] - m for m, c in zip(plm_minus_cons, plm_minus_cons_ci)]

colors_c = ['#2E7D32', '#2E7D32', '#757575', '#2E7D32', '#757575', '#757575']
edge_c = ['#1B5E20', '#1B5E20', '#424242', '#1B5E20', '#424242', '#424242']

x_c = np.arange(len(contrast_pairs))
ax2.bar(x_c, plm_minus_cons, width=0.48, yerr=[err_lo_c, err_hi_c], capsize=4,
        color=colors_c, edgecolor=edge_c, alpha=0.88)
ax2.axhline(0.0, color='black', linestyle='-', linewidth=0.9)

ax2.set_ylabel('Paired Block AUROC Difference: PLM - Conservation', fontsize=10.5, fontweight='bold')
ax2.set_title('B. Exploratory Out-of-Block Generalization Advantage', fontsize=11.5, fontweight='bold', pad=10)
ax2.set_xticks(x_c)
ax2.set_xticklabels(contrast_pairs, fontsize=8.2, rotation=20, ha='right')
ax2.set_ylim(-0.12, 0.35)
ax2.grid(True, axis='y', alpha=0.6)

# Labels for CI excludes zero vs spans zero
for i in range(len(contrast_pairs)):
    excludes = plm_minus_cons_ci[i][0] > 0
    stat_str = "Excludes 0" if excludes else "Spans 0"
    ax2.text(x_c[i], max(plm_minus_cons[i], 0) + err_hi_c[i] + 0.02, f"{plm_minus_cons[i]:+.3f}\n({stat_str})",
             ha='center', fontsize=7.2, fontweight='bold', color='#1B5E20' if excludes else '#424242')

ax2.annotate('Target dissociation:\nPLMs superior on allosteric UL30-ACV\nIndistinguishable on catalytic UL23',
             xy=(0, plm_minus_cons[0]), xytext=(0.8, 0.26),
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.1),
             fontsize=8.0, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#F1F8E9", ec="#33691E"))

plt.tight_layout()
save_fig(fig, "figure4_feature_family_degradation_matrix.png", ["figure4_temporal_shift_calibration.png"])

print("All publication figures 1-4 successfully generated!")
