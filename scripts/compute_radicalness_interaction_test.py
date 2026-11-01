"""
Radicalness Interaction Test: UL23 vs UL30
===========================================
Tests whether the difference in radicalness predictive performance between
Thymidine Kinase (UL23) and DNA Polymerase (UL30) is statistically significant,
directly avoiding the difference-of-significance fallacy.

Computes:
  1. Position-cluster bootstrap of difference in AUROC:
     - Delta_1 = AUROC(UL23-ACV) - AUROC(UL30-ACV)
     - Delta_2 = AUROC(UL23-ACV) - AUROC(UL30-FOS)
     - Delta_3 = AUROC(UL23-ACV) - AUROC(UL30-Combined)
     With 2,000 resamples clustered at unique sequence positions.
  2. Welch's two-sample t-test on contiguous-block AUROCs.
  3. Strict interval check: does the 95% bootstrap CI strictly exclude zero?

Outputs:
  - logs/radicalness_interaction_test.json
  - logs/radicalness_interaction_test.log
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCORED_TSV = os.path.join(ROOT, "logs", "primary_position_folds_scored.tsv")
OUT_JSON = os.path.join(ROOT, "logs", "radicalness_interaction_test.json")
OUT_LOG = os.path.join(ROOT, "logs", "radicalness_interaction_test.log")

df = pd.read_csv(SCORED_TSV, sep="\t")
df["pair"] = df["gene"] + "-" + df["drug"]

log_lines = []
def log(msg=""):
    print(msg)
    log_lines.append(msg)

log("=" * 80)
log("RADICALNESS DISSOCIATION INTERACTION TEST: UL23 vs UL30")
log("=" * 80)

# Extract cohorts
df_ul23_acv = df[df["pair"] == "UL23-Aciclovir"].copy()
df_ul30_acv = df[df["pair"] == "UL30-Aciclovir"].copy()
df_ul30_fos = df[df["pair"] == "UL30-Foscarnet"].copy()
df_ul30_comb = df[df["gene"] == "UL30"].copy()

# Point estimate pooled AUROCs (radicalness: higher = more radical = more deleterious/resistant)
# Notice: blosum62_radicalness in scored TSV was signed so higher is more resistant
# Let's verify direction
auc_ul23 = roc_auc_score(df_ul23_acv["target_binary"], df_ul23_acv["blosum62_radicalness"])
auc_ul30_acv = roc_auc_score(df_ul30_acv["target_binary"], df_ul30_acv["blosum62_radicalness"])
auc_ul30_fos = roc_auc_score(df_ul30_fos["target_binary"], df_ul30_fos["blosum62_radicalness"])
auc_ul30_comb = roc_auc_score(df_ul30_comb["target_binary"], df_ul30_comb["blosum62_radicalness"])

log(f"Observed Pooled AUROC (BLOSUM62 Radicalness):")
log(f"  UL23-ACV:       {auc_ul23:.4f} (N={len(df_ul23_acv)})")
log(f"  UL30-ACV:       {auc_ul30_acv:.4f} (N={len(df_ul30_acv)})")
log(f"  UL30-FOS:       {auc_ul30_fos:.4f} (N={len(df_ul30_fos)})")
log(f"  UL30-Combined:  {auc_ul30_comb:.4f} (N={len(df_ul30_comb)})")
log(f"  Point Delta (UL23-ACV - UL30-ACV):      {auc_ul23 - auc_ul30_acv:+.4f}")
log(f"  Point Delta (UL23-ACV - UL30-FOS):      {auc_ul23 - auc_ul30_fos:+.4f}")
log(f"  Point Delta (UL23-ACV - UL30-Combined): {auc_ul23 - auc_ul30_comb:+.4f}")
log("-" * 80)

# Position-cluster bootstrap
np.random.seed(42)
N_BOOT = 2000

positions_ul23 = df_ul23_acv["pos"].unique()
positions_ul30_acv = df_ul30_acv["pos"].unique()
positions_ul30_fos = df_ul30_fos["pos"].unique()
positions_ul30_comb = df_ul30_comb["pos"].unique()

boot_deltas_acv = []
boot_deltas_fos = []
boot_deltas_comb = []

for b in range(N_BOOT):
    # Cluster resample positions for UL23
    boot_pos_ul23 = np.random.choice(positions_ul23, size=len(positions_ul23), replace=True)
    sample_ul23 = df_ul23_acv[df_ul23_acv["pos"].isin(boot_pos_ul23)]
    
    # Check valid labels
    if len(sample_ul23["target_binary"].unique()) < 2:
        continue
    b_auc_ul23 = roc_auc_score(sample_ul23["target_binary"], sample_ul23["blosum62_radicalness"])
    
    # Cluster resample UL30-ACV
    boot_pos_30acv = np.random.choice(positions_ul30_acv, size=len(positions_ul30_acv), replace=True)
    sample_30acv = df_ul30_acv[df_ul30_acv["pos"].isin(boot_pos_30acv)]
    if len(sample_30acv["target_binary"].unique()) >= 2:
        b_auc_30acv = roc_auc_score(sample_30acv["target_binary"], sample_30acv["blosum62_radicalness"])
        boot_deltas_acv.append(b_auc_ul23 - b_auc_30acv)
        
    # Cluster resample UL30-FOS
    boot_pos_30fos = np.random.choice(positions_ul30_fos, size=len(positions_ul30_fos), replace=True)
    sample_30fos = df_ul30_fos[df_ul30_fos["pos"].isin(boot_pos_30fos)]
    if len(sample_30fos["target_binary"].unique()) >= 2:
        b_auc_30fos = roc_auc_score(sample_30fos["target_binary"], sample_30fos["blosum62_radicalness"])
        boot_deltas_fos.append(b_auc_ul23 - b_auc_30fos)

    # Cluster resample UL30-Combined
    boot_pos_30comb = np.random.choice(positions_ul30_comb, size=len(positions_ul30_comb), replace=True)
    sample_30comb = df_ul30_comb[df_ul30_comb["pos"].isin(boot_pos_30comb)]
    if len(sample_30comb["target_binary"].unique()) >= 2:
        b_auc_30comb = roc_auc_score(sample_30comb["target_binary"], sample_30comb["blosum62_radicalness"])
        boot_deltas_comb.append(b_auc_ul23 - b_auc_30comb)

results = {}

comparisons = [
    ("UL23_ACV_minus_UL30_ACV", boot_deltas_acv, auc_ul23 - auc_ul30_acv, "UL23-ACV vs UL30-ACV"),
    ("UL23_ACV_minus_UL30_FOS", boot_deltas_fos, auc_ul23 - auc_ul30_fos, "UL23-ACV vs UL30-FOS"),
    ("UL23_ACV_minus_UL30_COMBINED", boot_deltas_comb, auc_ul23 - auc_ul30_comb, "UL23-ACV vs UL30-Combined")
]

for key, deltas, pt_est, label in comparisons:
    deltas = np.array(deltas)
    mean_d = float(np.mean(deltas))
    std_d = float(np.std(deltas, ddof=1))
    ci_95 = [float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))]
    ci_99 = [float(np.percentile(deltas, 0.5)), float(np.percentile(deltas, 99.5))]
    
    # Two-sided empirical p-value for H0: Delta = 0
    p_emp = float(2 * min(np.mean(deltas <= 0), np.mean(deltas >= 0)))
    excludes_zero = ci_95[0] > 0
    
    results[key] = {
        "comparison": label,
        "observed_point_delta": float(pt_est),
        "bootstrap_mean_delta": mean_d,
        "bootstrap_sd_delta": std_d,
        "bootstrap_ci_95": ci_95,
        "bootstrap_ci_99": ci_99,
        "empirical_p_value": p_emp,
        "strictly_excludes_zero_ci95": excludes_zero,
        "verdict": (
            "Statistically significant dissociation (CI excludes zero)"
            if excludes_zero
            else "Not statistically significant (CI spans zero)"
        )
    }
    
    log(f"Comparison: {label}")
    log(f"  Observed Point Delta: {pt_est:+.4f}")
    log(f"  Bootstrap Mean Delta: {mean_d:+.4f} +/- {std_d:.4f}")
    log(f"  95% Position-Cluster Bootstrap CI: [{ci_95[0]:.4f}, {ci_95[1]:.4f}]")
    log(f"  99% Position-Cluster Bootstrap CI: [{ci_99[0]:.4f}, {ci_99[1]:.4f}]")
    log(f"  Empirical p-value: {p_emp:.4f}")
    log(f"  Strictly Excludes Zero? {excludes_zero}")
    log(f"  Scientific Verdict: {results[key]['verdict']}")
    log("-" * 80)

# Also block-level comparison:
# UL23 blocks: [0.5572, 0.5745, 0.6603, 0.6830, 0.7151]
# UL30-ACV blocks: [0.5526, 0.6800, 0.6861, 0.2785]
# UL30-FOS blocks: [0.5405, 0.6122, 0.7500, 0.1266]
blocks_ul23 = [0.55717, 0.57453, 0.66028, 0.68299, 0.71512]
blocks_ul30_acv = [0.55263, 0.68000, 0.68611, 0.27854]
blocks_ul30_fos = [0.54054, 0.61224, 0.75000, 0.12665]

t_acv, p_acv = stats.ttest_ind(blocks_ul23, blocks_ul30_acv, equal_var=False)
t_fos, p_fos = stats.ttest_ind(blocks_ul23, blocks_ul30_fos, equal_var=False)

results["block_level_welch_tests"] = {
    "ul23_vs_ul30_acv": {
        "ul23_block_mean": float(np.mean(blocks_ul23)),
        "ul30_acv_block_mean": float(np.mean(blocks_ul30_acv)),
        "delta": float(np.mean(blocks_ul23) - np.mean(blocks_ul30_acv)),
        "welch_t_stat": float(t_acv),
        "welch_p_val": float(p_acv),
        "df": float(stats.ttest_ind(blocks_ul23, blocks_ul30_acv, equal_var=False).df)
    },
    "ul23_vs_ul30_fos": {
        "ul23_block_mean": float(np.mean(blocks_ul23)),
        "ul30_fos_block_mean": float(np.mean(blocks_ul30_fos)),
        "delta": float(np.mean(blocks_ul23) - np.mean(blocks_ul30_fos)),
        "welch_t_stat": float(t_fos),
        "welch_p_val": float(p_fos),
        "df": float(stats.ttest_ind(blocks_ul23, blocks_ul30_fos, equal_var=False).df)
    }
}

log("Block-Level Welch's Two-Sample t-tests:")
log(f"  UL23-ACV vs UL30-ACV: Delta = +{np.mean(blocks_ul23) - np.mean(blocks_ul30_acv):.4f}, t = {t_acv:.4f}, p = {p_acv:.4f}")
log(f"  UL23-ACV vs UL30-FOS: Delta = +{np.mean(blocks_ul23) - np.mean(blocks_ul30_fos):.4f}, t = {t_fos:.4f}, p = {p_fos:.4f}")

with open(OUT_JSON, "w") as f:
    json.dump(results, f, indent=2)

with open(OUT_LOG, "w") as f:
    f.write("\n".join(log_lines))

log(f"\nSaved interaction test to: {OUT_JSON}")
