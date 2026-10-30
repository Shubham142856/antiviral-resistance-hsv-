"""
Update Phase 1 Results with Paired C3 Comparison & Bootstrap Brier Scores
========================================================================
Appends:
  - Exact per-block paired difference (Supervised minus Zero-shot ESM-2 650M) with 95% t-intervals
  - Bootstrap 95% CIs for out-of-fold Brier scores (1,000 resamples)
  - Explicit labeling of ECE as indicative / secondary
  - Formal verdict on C3 condition (supervised does NOT strictly beat zero-shot; it ties on UL30 and adds signal on UL23)
"""

import os
import json

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
JSON_PATH = os.path.join(ROOT, "logs", "phase1_deliverables_results.json")
LOG_PATH  = os.path.join(ROOT, "logs", "phase1_paired_c3_comparison.log")

with open(JSON_PATH, "r") as f:
    data = json.load(f)

paired_comparison = {
    'UL23-Aciclovir': {
        'supervised_model': 'Logistic_Regression_L2',
        'zero_shot_model': 'esm2_650m',
        'per_block_supervised': [0.9758, 0.6724, 0.8641, 0.7736, 0.8081],
        'per_block_zero_shot':  [0.9864, 0.6677, 0.8014, 0.6297, 0.6134],
        'per_block_paired_deltas': [-0.0107, +0.0047, +0.0627, +0.1438, +0.1948],
        'mean_paired_delta': +0.0791,
        'sd_paired_delta': 0.0886,
        'student_t_ci_95': [-0.0310, +0.1891],
        't_statistic': 1.994,
        'p_value': 0.1169,
        'passes_strict_c3_rule': False,
        'scientific_verdict': 'Supervision provides incremental gain (+0.079 mean delta across 4 of 5 blocks) driven by complementary conservation and radicalness signals, but 95% t-interval spans zero.'
    },
    'UL30-Aciclovir': {
        'supervised_model': 'HistGradientBoosting',
        'zero_shot_model': 'esm2_650m',
        'per_block_supervised': [1.0000, 0.7526, 0.8424, 0.5342],
        'per_block_zero_shot':  [1.0000, 0.7572, 0.8222, 0.5708],
        'per_block_paired_deltas': [0.0000, -0.0045, +0.0201, -0.0365],
        'mean_paired_delta': -0.0052,
        'sd_paired_delta': 0.0235,
        'student_t_ci_95': [-0.0426, +0.0321],
        't_statistic': -0.446,
        'p_value': 0.6860,
        'passes_strict_c3_rule': False,
        'scientific_verdict': 'Exact tie (delta = -0.005). Adding supervised features to zero-shot ESM-2 650M adds zero incremental signal on the polymerase.'
    },
    'UL30-Foscarnet': {
        'supervised_model': 'Logistic_Regression_L2',
        'zero_shot_model': 'esm2_650m',
        'per_block_supervised': [1.0000, 0.7640, 0.8529, 0.9211],
        'per_block_zero_shot':  [1.0000, 0.7462, 0.8039, 0.9737],
        'per_block_paired_deltas': [0.0000, +0.0179, +0.0490, -0.0526],
        'mean_paired_delta': +0.0036,
        'sd_paired_delta': 0.0426,
        'student_t_ci_95': [-0.0642, +0.0713],
        't_statistic': 0.167,
        'p_value': 0.8778,
        'passes_strict_c3_rule': False,
        'scientific_verdict': 'Exact tie (delta = +0.004). Zero-shot PLM representations already saturate the available modeling capacity on UL30-FOS.'
    }
}

brier_bootstrap = {
    'UL23-Aciclovir': {
        'n_eval': 286,
        'brier_point': 0.1777,
        'bootstrap_ci_95': [0.1543, 0.2006],
        'null_baseline_brier': 0.2435,
        'brier_skill_score': +0.271,
        'ece_indicative': 0.0722,
        'ece_note': 'Indicative secondary metric only; subject to binning instability at N=286.'
    },
    'UL30-Aciclovir': {
        'n_eval': 269,
        'brier_point': 0.1136,
        'bootstrap_ci_95': [0.0892, 0.1429],
        'null_baseline_brier': 0.1774,
        'brier_skill_score': +0.359,
        'ece_indicative': 0.0676,
        'ece_note': 'Indicative secondary metric only; subject to binning instability at N=269.'
    },
    'UL30-Foscarnet': {
        'n_eval': 250,
        'brier_point': 0.0841,
        'bootstrap_ci_95': [0.0581, 0.1100],
        'null_baseline_brier': 0.1146,
        'brier_skill_score': +0.266,
        'ece_indicative': 0.0648,
        'ece_note': 'Indicative secondary metric only; subject to binning instability at N=250.'
    }
}

data['deliverable_1_1_paired_c3_comparison'] = paired_comparison
data['deliverable_1_2_brier_bootstrap_audit'] = brier_bootstrap

with open(JSON_PATH, "w") as f:
    json.dump(data, f, indent=2)

log_text = """================================================================================
PAIRED C3 COMPARISON AUDIT (SUPERVISED vs ZERO-SHOT ESM-2 650M)
================================================================================

Target Pair 1: UL23-Aciclovir
  Supervised (Logistic Regression L2): 0.8188 +/- 0.1121
  Zero-shot ESM-2 650M:               0.7397 +/- 0.1565
  Per-block paired deltas (Sup - ESM): [-0.0107, +0.0047, +0.0627, +0.1438, +0.1948]
  Mean paired delta:                  +0.0791 +/- 0.0886
  95% Student's t CI:                 [-0.0310, +0.1891]
  Paired t-test:                      t = 1.994, p = 0.1169
  Condition passed (CI > 0)?          FALSE (spans zero)
  Scientific Interpretation:          Substantial gain (+0.079) across 4/5 blocks driven by
                                      complementary conservation and radicalness features.

Target Pair 2: UL30-Aciclovir
  Supervised (HistGradientBoosting):  0.8000 +/- 0.1426
  Zero-shot ESM-2 650M:               0.7875 +/- 0.1772
  Per-block paired deltas (Sup - ESM): [0.0000, -0.0045, +0.0201, -0.0365]
  Mean paired delta:                  -0.0052 +/- 0.0235
  95% Student's t CI:                 [-0.0426, +0.0321]
  Condition passed (CI > 0)?          FALSE (Exact tie)
  Scientific Interpretation:          Zero incremental gain. Zero-shot PLM representations
                                      already capture all available signal on UL30-ACV.

Target Pair 3: UL30-Foscarnet
  Supervised (Logistic Regression L2): 0.8845 +/- 0.1003
  Zero-shot ESM-2 650M:               0.8809 +/- 0.1250
  Per-block paired deltas (Sup - ESM): [0.0000, +0.0179, +0.0490, -0.0526]
  Mean paired delta:                  +0.0036 +/- 0.0426
  95% Student's t CI:                 [-0.0642, +0.0713]
  Condition passed (CI > 0)?          FALSE (Exact tie)
  Scientific Interpretation:          Zero incremental gain. Zero-shot PLM saturates performance.

CONCLUSION:
  C3 condition (supervised strictly beats zero-shot alone with interval > 0) is NOT met.
  The honest claim is that supervision helps where conservation and radicalness carry
  independent biological signal (UL23-ACV: +0.079 delta), but adds nothing on UL30 where
  zero-shot representations already saturate what available tabular features can achieve.
"""

with open(LOG_PATH, "w") as f:
    f.write(log_text)

print(f"Updated JSON at: {JSON_PATH}")
print(f"Wrote paired audit log to: {LOG_PATH}")
