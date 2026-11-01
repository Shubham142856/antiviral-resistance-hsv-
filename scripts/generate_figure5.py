"""
Generate Publication Figure 5: 3D Spatial-Neighborhood Holdout & Structural Generalization
========================================================================================
Visualizes:
  - Panel A: 3D C-alpha coordinate clustering into 5 spatial neighborhoods (KMeans on cryo-EM 8V1Q / crystal 1OF1).
  - Panel B: Dual-Regime Model Generalization (1D Contiguous-Block vs 3D Spatial-Neighborhood Holdout).
  - Panel C: Biological Target Dissociation Invariance across 1D and 3D Spatial Regimes.

Saves 300-DPI publication PNGs to:
  1. manuscript/figures/figure5_spatial_3d_generalization.png
  2. C:/Users/HP/.gemini/antigravity/brain/7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17/figures/figure5_spatial_3d_generalization.png
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Publication aesthetic styling
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['grid.color'] = '#EAEAEA'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.linewidth'] = 0.5

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SPATIAL_TSV = os.path.join(ROOT, "data", "processed", "primary_3d_spatial_folds.tsv")
RESULTS_JSON = os.path.join(ROOT, "logs", "spatial_3d_holdout_results.json")

ARTIFACT_DIR = r"C:\Users\HP\.gemini\antigravity\brain\7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17\figures"
WORKSPACE_DIR = os.path.join(ROOT, "manuscript", "figures")

os.makedirs(ARTIFACT_DIR, exist_ok=True)
os.makedirs(WORKSPACE_DIR, exist_ok=True)

df = pd.read_csv(SPATIAL_TSV, sep='\t')
with open(RESULTS_JSON, 'r') as f:
    res = json.load(f)

fig = plt.figure(figsize=(15, 5.2), dpi=300)

# -----------------------------------------------------------------------------
# Panel A: 3D Cryo-EM Spatial Clusters (UL30, 8V1Q)
# -----------------------------------------------------------------------------
ax1 = fig.add_subplot(1, 3, 1, projection='3d')
ul30_sub = df[(df['gene'] == 'UL30') & (df['has_3d_coord'] == 1)].drop_duplicates(subset=['pos'])

cluster_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
for c_id in range(5):
    c_pts = ul30_sub[ul30_sub['spatial_cluster_id'] == c_id]
    ax1.scatter(
        c_pts['coord_x'], c_pts['coord_y'], c_pts['coord_z'],
        c=cluster_colors[c_id], s=18, alpha=0.75,
        label=f"Fold {c_id} (N={len(c_pts)})", edgecolors='none'
    )

ax1.set_title("A  Cryo-EM Spatial Partitioning (UL30, 8V1Q)", fontsize=11, fontweight='bold', pad=8)
ax1.set_xlabel("X (Å)", fontsize=8, labelpad=-4)
ax1.set_ylabel("Y (Å)", fontsize=8, labelpad=-4)
ax1.set_zlabel("Z (Å)", fontsize=8, labelpad=-4)
ax1.tick_params(labelsize=7, pad=-2)
ax1.view_init(elev=22, azim=45)
ax1.legend(loc='upper right', fontsize=7, framealpha=0.8, handletextpad=0.2, borderaxespad=0.2)

# -----------------------------------------------------------------------------
# Panel B: 1D Block vs 3D Spatial AUROC (UL30-ACV & UL30-FOS)
# -----------------------------------------------------------------------------
ax2 = fig.add_subplot(1, 3, 2)

models = ['ESM-2 650M', 'ProtBERT', 'PSSM', 'Entropy', 'BLOSUM62']
m_keys = [
    'ESM-2 650M (zero-shot)',
    'ProtBERT (zero-shot)',
    'Multi-Herpesvirus PSSM',
    'Positional Shannon Entropy (-Entropy)',
    'BLOSUM62 Radicalness'
]

# Extract data for UL30-ACV
data_acv = res['results_by_pair']['UL30-Aciclovir']['resolved_core']
auc_1d = [data_acv[k]['1d_block']['mean_auc'] for k in m_keys]
err_1d = [data_acv[k]['1d_block']['sem_auc'] for k in m_keys]
auc_3d = [data_acv[k]['3d_spatial']['mean_auc'] for k in m_keys]
err_3d = [data_acv[k]['3d_spatial']['sem_auc'] for k in m_keys]

x = np.arange(len(models))
width = 0.35

rects1 = ax2.bar(x - width/2, auc_1d, width, yerr=err_1d, capsize=3,
                 color='#4A90E2', edgecolor='#1F4E79', linewidth=0.8, label='1D Contiguous Block')
rects2 = ax2.bar(x + width/2, auc_3d, width, yerr=err_3d, capsize=3,
                 color='#50E3C2', edgecolor='#1B7A66', linewidth=0.8, label='3D Spatial Neighborhood')

ax2.axhline(0.50, color='#999999', linestyle=':', linewidth=1.0, label='Chance (0.50)')
ax2.set_title("B  Cross-Regime Stability (UL30-ACV)", fontsize=11, fontweight='bold', pad=10)
ax2.set_ylabel("Mean Cross-Validation AUROC", fontsize=9)
ax2.set_xticks(x)
ax2.set_xticklabels(models, rotation=25, ha='right', fontsize=8)
ax2.set_ylim(0.40, 0.95)
ax2.grid(True, axis='y', alpha=0.6)
ax2.legend(loc='upper right', fontsize=8, framealpha=0.9)

# -----------------------------------------------------------------------------
# Panel C: Radicalness Dissociation Invariance (UL23 vs UL30)
# -----------------------------------------------------------------------------
ax3 = fig.add_subplot(1, 3, 3)

pairs_labels = ['UL23-ACV\n(TK, Kinase)', 'UL30-ACV\n(Pol, Allosteric)', 'UL30-FOS\n(Pol, Allosteric)']
pairs_keys = ['UL23-Aciclovir', 'UL30-Aciclovir', 'UL30-Foscarnet']

rad_1d = [res['results_by_pair'][pk]['resolved_core']['BLOSUM62 Radicalness']['1d_block']['mean_auc'] for pk in pairs_keys]
rad_1d_err = [res['results_by_pair'][pk]['resolved_core']['BLOSUM62 Radicalness']['1d_block']['sem_auc'] for pk in pairs_keys]

rad_3d = [res['results_by_pair'][pk]['resolved_core']['BLOSUM62 Radicalness']['3d_spatial']['mean_auc'] for pk in pairs_keys]
rad_3d_err = [res['results_by_pair'][pk]['resolved_core']['BLOSUM62 Radicalness']['3d_spatial']['sem_auc'] for pk in pairs_keys]

x3 = np.arange(len(pairs_labels))
rects3_1 = ax3.bar(x3 - width/2, rad_1d, width, yerr=rad_1d_err, capsize=3,
                   color='#E74C3C', edgecolor='#922B21', linewidth=0.8, label='BLOSUM62 (1D Block)')
rects3_2 = ax3.bar(x3 + width/2, rad_3d, width, yerr=rad_3d_err, capsize=3,
                   color='#F39C12', edgecolor='#B9770E', linewidth=0.8, label='BLOSUM62 (3D Spatial)')

ax3.axhline(0.50, color='#999999', linestyle=':', linewidth=1.0, label='Chance (0.50)')
ax3.set_title("C  Target Dissociation Invariance", fontsize=11, fontweight='bold', pad=10)
ax3.set_ylabel("BLOSUM62 Radicalness AUROC", fontsize=9)
ax3.set_xticks(x3)
ax3.set_xticklabels(pairs_labels, fontsize=8)
ax3.set_ylim(0.40, 0.75)
ax3.grid(True, axis='y', alpha=0.6)
ax3.legend(loc='upper right', fontsize=8, framealpha=0.9)

# Value annotations on bars
for bar in rects3_1:
    h = bar.get_height()
    ax3.text(bar.get_x() + bar.get_width()/2., h + 0.015, f"{h:.3f}", ha='center', va='bottom', fontsize=7, fontweight='bold')
for bar in rects3_2:
    h = bar.get_height()
    ax3.text(bar.get_x() + bar.get_width()/2., h + 0.015, f"{h:.3f}", ha='center', va='bottom', fontsize=7, fontweight='bold')

plt.tight_layout()

# Save
out_art = os.path.join(ARTIFACT_DIR, "figure5_spatial_3d_generalization.png")
out_ws  = os.path.join(WORKSPACE_DIR, "figure5_spatial_3d_generalization.png")

plt.savefig(out_art, dpi=300, bbox_inches='tight')
plt.savefig(out_ws, dpi=300, bbox_inches='tight')
plt.close()

print(f"Successfully generated Figure 5:")
print(f"  Artifact:  {out_art}")
print(f"  Workspace: {out_ws}")

