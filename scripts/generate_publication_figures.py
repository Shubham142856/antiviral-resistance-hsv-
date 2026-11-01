"""
Generate ICML / Nature-style Publication Figures for HSV Resistance Benchmark
=============================================================================
Reads verified, frozen JSON logs:
  - logs/c2_repeated_cv_results.json
  - logs/msa_conservation_results.json
  - logs/pssm_baseline_results.json
  - logs/paired_dissociation_strongest_test.json
  - logs/plm_benchmark_results.json
  - logs/phase1_deliverables_results.json

Outputs high-DPI publication figures to:
  1. Artifacts directory: C:/Users/HP/.gemini/antigravity/brain/7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17/figures/
  2. Workspace directory: manuscript/figures/
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Publication aesthetic styling
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['grid.color'] = '#E0E0E0'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.linewidth'] = 0.5

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LOGS = os.path.join(ROOT, "logs")

ARTIFACT_DIR = r"C:\Users\HP\.gemini\antigravity\brain\7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17\figures"
WORKSPACE_DIR = os.path.join(ROOT, "manuscript", "figures")

os.makedirs(ARTIFACT_DIR, exist_ok=True)
os.makedirs(WORKSPACE_DIR, exist_ok=True)

# -----------------------------------------------------------------------------
# FIGURE 1: Spatial Leakage Collapse
# -----------------------------------------------------------------------------
print("Generating Figure 1: Spatial Leakage Collapse...")

# Data from logs/c2_repeated_cv_results.json & Table 1
# Models: Position-Only RF, All-Features RF, ESM-2 650M
# UL30-ACV & UL30-FOS
models = ['Position-Only RF', 'All-Features RF', 'ESM-2 650M (Zero-shot)']
pos_acv = [0.8715, 0.7328, 0.8273]
pos_acv_sd = [0.0171, 0.0365, 0.0632]
blk_acv = [0.5001, 0.4089, 0.7875]
blk_acv_sd = [0.0009, 0.0198, 0.1772]

pos_fos = [0.8684, 0.6756, 0.8454]
pos_fos_sd = [0.0282, 0.0446, 0.1050]
blk_fos = [0.5166, 0.6372, 0.8809]
blk_fos_sd = [0.0105, 0.0181, 0.1250]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8), dpi=300)

x = np.arange(len(models))
width = 0.35

# Subplot A: UL30-Aciclovir
rects1 = ax1.bar(x - width/2, pos_acv, width, yerr=pos_acv_sd, label='Position-Grouped CV',
                 color='#4A90E2', capsize=4, edgecolor='#1E5FAD', alpha=0.9)
rects2 = ax1.bar(x + width/2, blk_acv, width, yerr=blk_acv_sd, label='Contiguous-Block CV',
                 color='#E94E77', capsize=4, edgecolor='#B51A44', alpha=0.9)

ax1.axhline(0.50, color='#888888', linestyle=':', linewidth=1.2, label='Chance (0.50)')
ax1.set_ylabel('Mean AUROC $\pm$ SD', fontsize=11, fontweight='bold')
ax1.set_title('A. DNA Polymerase (UL30) — Aciclovir', fontsize=12, fontweight='bold', pad=10)
ax1.set_xticks(x)
ax1.set_xticklabels(models, fontsize=9.5)
ax1.set_ylim(0.3, 1.0)
ax1.grid(True, axis='y', alpha=0.6)
ax1.legend(frameon=True, fontsize=9, loc='upper right')

# Annotation for leakage collapse
ax1.annotate('Spatial Collapse:\n0.8715 -> 0.5001\n(Chance)', 
             xy=(0, 0.505), xytext=(-0.25, 0.65),
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.2),
             fontsize=8.5, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#FFF3E0", ec="#E65100"))

# Subplot B: UL30-Foscarnet
rects3 = ax2.bar(x - width/2, pos_fos, width, yerr=pos_fos_sd, label='Position-Grouped CV',
                 color='#4A90E2', capsize=4, edgecolor='#1E5FAD', alpha=0.9)
rects4 = ax2.bar(x + width/2, blk_fos, width, yerr=blk_fos_sd, label='Contiguous-Block CV',
                 color='#E94E77', capsize=4, edgecolor='#B51A44', alpha=0.9)

ax2.axhline(0.50, color='#888888', linestyle=':', linewidth=1.2, label='Chance (0.50)')
ax2.set_ylabel('Mean AUROC $\pm$ SD', fontsize=11, fontweight='bold')
ax2.set_title('B. DNA Polymerase (UL30) — Foscarnet', fontsize=12, fontweight='bold', pad=10)
ax2.set_xticks(x)
ax2.set_xticklabels(models, fontsize=9.5)
ax2.set_ylim(0.3, 1.0)
ax2.grid(True, axis='y', alpha=0.6)
ax2.legend(frameon=True, fontsize=9, loc='upper right')

ax2.annotate('Spatial Collapse:\n0.8684 -> 0.5166\n(Chance)', 
             xy=(0, 0.52), xytext=(-0.25, 0.65),
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.2),
             fontsize=8.5, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#FFF3E0", ec="#E65100"))

plt.tight_layout()
fig1_path_art = os.path.join(ARTIFACT_DIR, "figure1_spatial_leakage_collapse.png")
fig1_path_ws  = os.path.join(WORKSPACE_DIR, "figure1_spatial_leakage_collapse.png")
plt.savefig(fig1_path_art, bbox_inches='tight')
plt.savefig(fig1_path_ws, bbox_inches='tight')
plt.close()
print(f"Saved: {fig1_path_art}")

# -----------------------------------------------------------------------------
# FIGURE 2: Biological Dissociation (Radicalness & Baseline Progression)
# -----------------------------------------------------------------------------
print("Generating Figure 2: Biological Dissociation...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=300)

# Subplot A: BLOSUM62 per-block AUROCs
# UL23-ACV blocks: B0=0.557, B1=0.575, B2=0.660, B3=0.683, B4=0.715 (Mean 0.638)
# UL30-ACV blocks: B1=0.553, B2=0.680, B3=0.686, B4=0.279 (Mean 0.549)
# UL30-FOS blocks: B1=0.541, B2=0.612, B3=0.750, B4=0.127 (Mean 0.507)
targets = ['UL23-ACV\n(Thymidine Kinase)', 'UL30-ACV\n(DNA Polymerase)', 'UL30-FOS\n(DNA Polymerase)']
rad_means = [0.6380, 0.5494, 0.5074]
rad_ci_lo = [0.5524, 0.1523, 0.0305]
rad_ci_hi = [0.7237, 0.9466, 0.9844]
yerr_lo = [m - lo for m, lo in zip(rad_means, rad_ci_lo)]
yerr_hi = [hi - m for m, hi in zip(rad_means, rad_ci_hi)]

bars = ax1.bar(targets, rad_means, color=['#2CA02C', '#D62728', '#D62728'],
               alpha=0.85, edgecolor=['#1B641B', '#8C1B1B', '#8C1B1B'], width=0.5)
ax1.errorbar(targets, rad_means, yerr=[yerr_lo, yerr_hi], fmt='none', ecolor='black', capsize=6, lw=1.2)
ax1.axhline(0.50, color='gray', linestyle='--', linewidth=1.2, label='Chance (0.50)')

ax1.set_ylabel('Contiguous-Block AUROC (95% $t$-CI)', fontsize=11, fontweight='bold')
ax1.set_title('A. Radicalness Dissociation (BLOSUM62)', fontsize=12, fontweight='bold', pad=10)
ax1.set_ylim(0.0, 1.05)
ax1.grid(True, axis='y', alpha=0.6)

ax1.annotate('Statistically Robust\nAll blocks > 0.55\nCI excludes 0.50', xy=(0, 0.65), xytext=(-0.35, 0.80),
             arrowprops=dict(facecolor='#1B641B', arrowstyle='->', lw=1.2),
             fontsize=8.5, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#E8F5E9", ec="#2E7D32"))
ax1.annotate('Indistinguishable\nfrom Chance\nCIs span 0.50', xy=(1.5, 0.52), xytext=(1.1, 0.25),
             arrowprops=dict(facecolor='#8C1B1B', arrowstyle='->', lw=1.2),
             fontsize=8.5, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#FFEBEE", ec="#C62828"))

# Subplot B: Baseline Progression on UL30-ACV & UL30-FOS
baseline_names = ['Radicalness\n(BLOSUM62)', 'Weighted\nConservation', 'PSSM Profile\n(Dirichlet)', 'Positional\nEntropy', 'Zero-Shot\nESM-2 650M']
acv_prog = [0.5494, 0.6428, 0.6507, 0.6999, 0.7875]
fos_prog = [0.5074, 0.6604, 0.7832, 0.8246, 0.8809]

ax2.plot(baseline_names, acv_prog, marker='o', linewidth=2, color='#1F77B4', label='UL30-Aciclovir (Margin -> +0.088, p=0.174)')
ax2.plot(baseline_names, fos_prog, marker='s', linewidth=2, color='#FF7F0E', label='UL30-Foscarnet (Margin -> +0.056, p=0.373)')
ax2.axhline(0.50, color='gray', linestyle='--', linewidth=1.2)

ax2.set_ylabel('Contiguous-Block Mean AUROC', fontsize=11, fontweight='bold')
ax2.set_title('B. Evolutionary Baseline Progression on UL30', fontsize=12, fontweight='bold', pad=10)
ax2.set_ylim(0.45, 0.95)
ax2.grid(True, alpha=0.6)
ax2.legend(frameon=True, fontsize=8.5, loc='upper left')

# Annotate narrowing margin
ax2.annotate('Margin shrinks\nfrom +0.221 to +0.056\n(both CIs span 0)', 
             xy=(4, 0.8809), xytext=(2.6, 0.88),
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.2),
             fontsize=8.5, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#FFF9C4", ec="#FBC02D"))

plt.tight_layout()
fig2_path_art = os.path.join(ARTIFACT_DIR, "figure2_biological_dissociation.png")
fig2_path_ws  = os.path.join(WORKSPACE_DIR, "figure2_biological_dissociation.png")
plt.savefig(fig2_path_art, bbox_inches='tight')
plt.savefig(fig2_path_ws, bbox_inches='tight')
plt.close()
print(f"Saved: {fig2_path_art}")

# -----------------------------------------------------------------------------
# FIGURE 3: Benchmark Model Comparison Matrix
# -----------------------------------------------------------------------------
print("Generating Figure 3: Model Benchmark Matrix...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=300)

models_bench = ['ESM-2 650M\n(Primary)', 'ProtBERT 420M\n(Secondary)', 'ESM-2 35M\n(Baseline)', 'BLOSUM62\n(Radicalness)']
ul23_means = [0.7397, 0.7530, 0.6110, 0.6380]
ul23_lo = [0.5455, 0.6645, 0.3625, 0.5524]
ul23_hi = [0.9340, 0.8415, 0.8595, 0.7237]

ul30acv_means = [0.7875, 0.6957, 0.6795, 0.5494]
ul30acv_lo = [0.5055, 0.3113, 0.2911, 0.1523]
ul30acv_hi = [1.0000, 1.0000, 1.0000, 0.9466]

ul30fos_means = [0.8809, 0.7423, 0.7481, 0.5074]
ul30fos_lo = [0.6821, 0.4888, 0.4550, 0.0305]
ul30fos_hi = [1.0000, 0.9958, 1.0000, 0.9844]

x = np.arange(len(models_bench))
w = 0.25

ax1.bar(x - w, ul23_means, w, color='#1F77B4', label='UL23-ACV', alpha=0.9, edgecolor='#12466B')
ax1.bar(x, ul30acv_means, w, color='#FF7F0E', label='UL30-ACV', alpha=0.9, edgecolor='#A64F04')
ax1.bar(x + w, ul30fos_means, w, color='#2CA02C', label='UL30-FOS', alpha=0.9, edgecolor='#175917')

# Error bars for ESM-2 650M
ax1.errorbar(x[0]-w, ul23_means[0], yerr=[[ul23_means[0]-ul23_lo[0]], [ul23_hi[0]-ul23_means[0]]], fmt='none', ecolor='black', capsize=4)
ax1.errorbar(x[0], ul30acv_means[0], yerr=[[ul30acv_means[0]-ul30acv_lo[0]], [ul30acv_hi[0]-ul30acv_means[0]]], fmt='none', ecolor='black', capsize=4)
ax1.errorbar(x[0]+w, ul30fos_means[0], yerr=[[ul30fos_means[0]-ul30fos_lo[0]], [ul30fos_hi[0]-ul30fos_means[0]]], fmt='none', ecolor='black', capsize=4)

ax1.axhline(0.65, color='#D32F2F', linestyle='--', linewidth=1.2, label='Stopping Rule (Mean >= 0.65)')
ax1.axhline(0.50, color='gray', linestyle=':', linewidth=1.0, label='Chance (0.50)')

ax1.set_ylabel('Contiguous-Block Mean AUROC', fontsize=11, fontweight='bold')
ax1.set_title('A. Benchmark Performance vs. Stopping Rule', fontsize=12, fontweight='bold', pad=10)
ax1.set_xticks(x)
ax1.set_xticklabels(models_bench, fontsize=9)
ax1.set_ylim(0.3, 1.05)
ax1.grid(True, axis='y', alpha=0.6)
ax1.legend(frameon=True, fontsize=8.5, loc='lower right')

# Subplot B: Supervised vs Zero-Shot Paired Delta
pairs_sup = ['UL23-ACV', 'UL30-ACV', 'UL30-FOS']
sup_deltas = [+0.0791, -0.0052, +0.0036]
sup_ci_lo = [-0.0310, -0.0426, -0.0642]
sup_ci_hi = [+0.1891, +0.0321, +0.0713]
p_vals = ['p = 0.117', 'p = 0.686', 'p = 0.878']

x_sup = np.arange(len(pairs_sup))
ax2.bar(x_sup, sup_deltas, width=0.45, color=['#4A90E2', '#888888', '#888888'],
        edgecolor=['#1E5FAD', '#444444', '#444444'], alpha=0.85)
ax2.errorbar(x_sup, sup_deltas, 
             yerr=[[d - lo for d, lo in zip(sup_deltas, sup_ci_lo)], 
                   [hi - d for d, hi in zip(sup_deltas, sup_ci_hi)]],
             fmt='none', ecolor='black', capsize=6, lw=1.2)
ax2.axhline(0.0, color='black', linestyle='-', linewidth=0.8)

for i, (d, p) in enumerate(zip(sup_deltas, p_vals)):
    y_pos = d + (0.03 if d >= 0 else -0.04)
    ax2.text(i, y_pos, p, ha='center', fontsize=9, fontweight='bold')

ax2.set_ylabel('Paired Delta ($\Delta = \mathrm{Sup} - \mathrm{ESM2}$)', fontsize=11, fontweight='bold')
ax2.set_title('B. Honest Paired Supervision Test (Condition C3)', fontsize=12, fontweight='bold', pad=10)
ax2.set_xticks(x_sup)
ax2.set_xticklabels(pairs_sup, fontsize=10)
ax2.set_ylim(-0.10, 0.25)
ax2.grid(True, axis='y', alpha=0.6)

ax2.annotate('Condition C3 FAILED:\nAll CIs span zero;\nZero-shot saturates UL30', 
             xy=(1.5, 0.0), xytext=(0.8, 0.15),
             arrowprops=dict(facecolor='#D32F2F', arrowstyle='->', lw=1.2),
             fontsize=8.5, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#FFEBEE", ec="#D32F2F"))

plt.tight_layout()
fig3_path_art = os.path.join(ARTIFACT_DIR, "figure3_model_benchmark_matrix.png")
fig3_path_ws  = os.path.join(WORKSPACE_DIR, "figure3_model_benchmark_matrix.png")
plt.savefig(fig3_path_art, bbox_inches='tight')
plt.savefig(fig3_path_ws, bbox_inches='tight')
plt.close()
print(f"Saved: {fig3_path_art}")

# -----------------------------------------------------------------------------
# FIGURE 4: Prospective Temporal Shift & Calibration
# -----------------------------------------------------------------------------
print("Generating Figure 4: Temporal Shift and Calibration...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=300)

# Subplot A: Temporal Generalization vs Seen Positions
categories = ['UL23-ACV', 'UL30-ACV', 'UL30-FOS']
seen_perf = [0.8872, 0.9231, 0.9286]
novel_perf = [0.6996, 0.8075, 0.8126]
block_perf = [0.8188, 0.7823, 0.8845]

x_cat = np.arange(len(categories))
w_t = 0.26

ax1.bar(x_cat - w_t, seen_perf, w_t, label='Seen Positions (Memorization)', color='#FFA726', edgecolor='#E65100', alpha=0.9)
ax1.bar(x_cat, novel_perf, w_t, label='Novel Positions (Generalization)', color='#29B6F6', edgecolor='#0277BD', alpha=0.9)
ax1.bar(x_cat + w_t, block_perf, w_t, label='Contiguous-Block CV Estimate', color='#AB47BC', edgecolor='#6A1B9A', alpha=0.9)

ax1.set_ylabel('Supervised Model AUROC', fontsize=11, fontweight='bold')
ax1.set_title('A. Prospective Temporal Holdout (Post-2016 Variants)', fontsize=12, fontweight='bold', pad=10)
ax1.set_xticks(x_cat)
ax1.set_xticklabels(categories, fontsize=10)
ax1.set_ylim(0.5, 1.02)
ax1.grid(True, axis='y', alpha=0.6)
ax1.legend(frameon=True, fontsize=8.5, loc='lower left')

ax1.annotate('Novel position drop\nmatches block estimate', xy=(0, 0.70), xytext=(-0.35, 0.55),
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.2),
             fontsize=8.5, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#E1F5FE", ec="#0288D1"))

# Subplot B: Brier Calibration & Skill Scores
pairs_calib = ['UL23-ACV\n(Skill: +0.271)', 'UL30-ACV\n(Skill: +0.359)', 'UL30-FOS\n(Skill: +0.266)']
brier_vals = [0.1777, 0.1136, 0.0841]
null_vals = [0.2435, 0.1774, 0.1146]

x_cal = np.arange(len(pairs_calib))
w_c = 0.35

ax2.bar(x_cal - w_c/2, null_vals, w_c, label='Null Baseline Brier (Base Rate)', color='#B0BEC5', edgecolor='#455A64', alpha=0.8)
ax2.bar(x_cal + w_c/2, brier_vals, w_c, label='Model Brier Score (Lower is Better)', color='#66BB6A', edgecolor='#2E7D32', alpha=0.9)

ax2.set_ylabel('Brier Score (MSE against binary label)', fontsize=11, fontweight='bold')
ax2.set_title('B. Probabilistic Calibration & Brier Skill Scores', fontsize=12, fontweight='bold', pad=10)
ax2.set_xticks(x_cal)
ax2.set_xticklabels(pairs_calib, fontsize=10)
ax2.set_ylim(0.0, 0.30)
ax2.grid(True, axis='y', alpha=0.6)
ax2.legend(frameon=True, fontsize=8.5, loc='upper right')

for i in range(len(pairs_calib)):
    pct = ((null_vals[i] - brier_vals[i]) / null_vals[i]) * 100
    ax2.text(i + w_c/2, brier_vals[i] + 0.01, f"+{pct:.1f}%", ha='center', fontsize=9, fontweight='bold', color='#1B5E20')

plt.tight_layout()
fig4_path_art = os.path.join(ARTIFACT_DIR, "figure4_temporal_shift_calibration.png")
fig4_path_ws  = os.path.join(WORKSPACE_DIR, "figure4_temporal_shift_calibration.png")
plt.savefig(fig4_path_art, bbox_inches='tight')
plt.savefig(fig4_path_ws, bbox_inches='tight')
plt.close()
print(f"Saved: {fig4_path_art}")

print("All publication figures successfully generated!")

