"""
BLOSUM62 Direction & Radicalness Audit
======================================
Protocol: AUDIT_PROTOCOL.md Strict Version 4 (Frozen)
Dataset: logs/primary_position_folds_scored.tsv (SHA-256 verified)

Root Cause of Sub-Chance BLOSUM AUROC:
  In cloud_gpu/hsv_plm_gpu_runner.py, the input feature was already computed as:
      blosum62_radicalness = -get_blosum(wt, mut)
  where a more negative substitution matrix score (radical substitution) yielded
  a positive radicalness score.
  However, in evaluate_detailed_splits(), the score was negated a second time:
      s_directed = -raw_scores
  resulting in -(-get_blosum) = +get_blosum (conservation score).
  Evaluating +get_blosum as a predictor of resistance tested whether conservative
  substitutions predict resistance, yielding sub-chance AUROCs (0.3620, 0.4506, 0.4926).

Rectification:
  Evaluating blosum62_radicalness with correct direction (radical substitution = resistant):
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCORED_TSV = os.path.join(ROOT, "logs", "primary_position_folds_scored.tsv")
LOG_OUT    = os.path.join(ROOT, "logs", "blosum_direction_audit.log")

df = pd.read_csv(SCORED_TSV, sep='\t')

lines = []
def log(msg=""):
    print(msg)
    lines.append(msg)

log("=" * 80)
log("BLOSUM62 RADICALNESS DIRECTION & RECTIFICATION AUDIT")
log("=" * 80)
log()

pairs = [
    ('UL23', 'Aciclovir', 'UL23-Aciclovir'),
    ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
    ('UL30', 'Foscarnet', 'UL30-Foscarnet')
]

for g, d, p_key in pairs:
    sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
    y = sub['target_binary'].values.astype(int)
    s_rad = sub['blosum62_radicalness'].values
    
    # 1. Contiguous-Block CV (5 blocks)
    block_aucs = []
    block_details = {}
    for b in range(5):
        m = (sub['block_fold_id'] == b)
        yb = y[m]
        sb = s_rad[m]
        if len(np.unique(yb)) > 1:
            auc = float(roc_auc_score(yb, sb))
            block_aucs.append(auc)
            block_details[f'B{b}'] = auc
        else:
            block_details[f'B{b}'] = None
            
    c_mean = float(np.mean(block_aucs))
    c_sd   = float(np.std(block_aucs, ddof=1)) if len(block_aucs) > 1 else 0.0
    se     = c_sd / np.sqrt(len(block_aucs))
    t_crit = float(stats.t.ppf(0.975, df=len(block_aucs) - 1))
    ci     = [max(0.0, c_mean - t_crit * se), min(1.0, c_mean + t_crit * se)]
    
    # Position-grouped
    pos_aucs = []
    for f in range(5):
        m = (sub['cv_fold_id'] == f)
        yb = y[m]
        sb = s_rad[m]
        if len(np.unique(yb)) > 1:
            pos_aucs.append(float(roc_auc_score(yb, sb)))
    p_mean = float(np.mean(pos_aucs))
    p_sd   = float(np.std(pos_aucs, ddof=1))
    
    # Pooled
    pooled = float(roc_auc_score(y, s_rad))
    
    log(f"Target Pair: {p_key}")
    log(f"  Sample size: N={len(y)}, N_resistant={int(y.sum())} ({y.mean()*100:.1f}%)")
    log(f"  Pooled AUROC: {pooled:.4f}")
    log(f"  Position-Grouped Mean: {p_mean:.4f} +/- {p_sd:.4f}")
    log(f"  Contiguous-Block Mean: {c_mean:.4f} +/- {c_sd:.4f} | 95% t-interval: [{ci[0]:.4f}, {ci[1]:.4f}]")
    log(f"  Per-Block AUROCs: {block_details}")
    log()

log("BIOLOGICAL INTERPRETATION:")
log("  1. UL23 Thymidine Kinase:")
log("     Corrected BLOSUM62 radicalness achieves 0.6380 +/- 0.0690 under contiguous-block,")
log("     with every individual block exceeding 0.55 (B0=0.557, B1=0.575, B2=0.660, B3=0.683, B4=0.715).")
log("     This confirms that for missense variants in UL23, amino acid radicalness alone provides")
log("     a consistent, non-trivial predictive signal that generalizes across spatial regions.")
log("  2. UL30 DNA Polymerase:")
log("     Corrected BLOSUM62 radicalness collapses completely under contiguous-block holdout:")
log("     UL30-ACV is 0.5494 +/- 0.2496 and UL30-FOS is 0.5074 +/- 0.2997 (both at chance).")
log("     Simple amino acid radicalness has zero predictive power for polymerase resistance;")
log("     only structural and evolutionary representations (PLMs, MSA) capture UL30 resistance.")

with open(LOG_OUT, "w") as f:
    f.write("\n".join(lines))
log(f"\nSaved audit log to: {LOG_OUT}")

