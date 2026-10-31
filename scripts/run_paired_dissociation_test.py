"""
Paired Dissociation Test: ESM-2 650M vs Negative MSA Positional Entropy
========================================================================
Computes:
  - Per-block paired delta (ESM-2 650M minus Negative MSA Entropy)
  - Mean paired delta
  - Standard deviation of paired delta (ddof=1)
  - Student's t-statistic, two-tailed p-value, and exact 95% t-interval
  - Strict falsification check: does the 95% t-interval strictly exclude zero?
"""

import os
import json
import math
import numpy as np
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PLM_JSON = os.path.join(ROOT, "logs", "plm_benchmark_results.json")
MSA_JSON = os.path.join(ROOT, "logs", "msa_conservation_results.json")
OUT_JSON = os.path.join(ROOT, "logs", "paired_dissociation_entropy_test.json")
OUT_LOG  = os.path.join(ROOT, "logs", "paired_dissociation_entropy_test.log")

with open(PLM_JSON, "r") as f:
    plm_data = json.load(f)

with open(MSA_JSON, "r") as f:
    msa_data = json.load(f)

pairs = ['UL23-Aciclovir', 'UL30-Aciclovir', 'UL30-Foscarnet']
results = {}
log_lines = []

def log(msg=""):
    print(msg)
    log_lines.append(msg)

log("=" * 85)
log("PAIRED DISSOCIATION TEST: ZERO-SHOT ESM-2 650M vs POSITIONAL ENTROPY")
log("=" * 85)

for pair in pairs:
    esm_blocks = plm_data['esm2_650m']['metrics'][pair]['per_block_breakdown']
    msa_blocks = msa_data['results_by_pair'][pair]['Negative MSA Entropy (-Entropy)']['per_block_auroc']
    
    paired_deltas = []
    block_names = []
    esm_vals = []
    msa_vals = []
    
    for b_id in ['B0', 'B1', 'B2', 'B3', 'B4']:
        e_val = esm_blocks[b_id]['auroc']
        m_val = msa_blocks.get(b_id)
        
        if e_val is not None and m_val is not None:
            block_names.append(b_id)
            esm_vals.append(e_val)
            msa_vals.append(m_val)
            delta = e_val - m_val
            paired_deltas.append(delta)
            
    n_eval = len(paired_deltas)
    mean_delta = float(np.mean(paired_deltas))
    sd_delta = float(np.std(paired_deltas, ddof=1)) if n_eval > 1 else 0.0
    se = sd_delta / math.sqrt(n_eval) if n_eval > 0 else 0.0
    
    df_t = n_eval - 1
    t_crit = float(stats.t.ppf(0.975, df=df_t))
    t_stat = mean_delta / se if se > 0 else 0.0
    p_val = float(2 * (1 - stats.t.cdf(abs(t_stat), df=df_t)))
    
    ci_95 = [mean_delta - t_crit * se, mean_delta + t_crit * se]
    excludes_zero = (ci_95[0] > 0 and ci_95[1] > 0) or (ci_95[0] < 0 and ci_95[1] < 0)
    
    results[pair] = {
        'evaluable_blocks': block_names,
        'n_eval': n_eval,
        'esm2_650m_per_block': esm_vals,
        'entropy_per_block': msa_vals,
        'esm2_650m_mean': float(np.mean(esm_vals)),
        'entropy_mean': float(np.mean(msa_vals)),
        'paired_deltas': paired_deltas,
        'mean_paired_delta': mean_delta,
        'sd_paired_delta': sd_delta,
        'se_paired_delta': se,
        't_statistic': t_stat,
        'degrees_of_freedom': df_t,
        'p_value': p_val,
        'student_t_ci_95': ci_95,
        'strictly_excludes_zero': excludes_zero,
        'scientific_verdict': (
            'Statistically significant advantage' if excludes_zero
            else 'Qualitative / directional advantage only (95% CI spans zero)'
        )
    }
    
    log(f"Target Pair: {pair}")
    log(f"  Evaluable Blocks: {block_names} (N={n_eval})")
    log(f"  ESM-2 650M per-block: {[round(x, 4) for x in esm_vals]} (Mean = {np.mean(esm_vals):.4f})")
    log(f"  Entropy    per-block: {[round(x, 4) for x in msa_vals]} (Mean = {np.mean(msa_vals):.4f})")
    log(f"  Paired deltas (ESM - Entropy): {[round(x, 4) for x in paired_deltas]}")
    log(f"  Mean paired delta: {mean_delta:+.4f} +/- {sd_delta:.4f} (SE = {se:.4f})")
    log(f"  t({df_t}) = {t_stat:.4f}, p = {p_val:.4f}")
    log(f"  95% Student's t CI: [{ci_95[0]:.4f}, {ci_95[1]:.4f}]")
    log(f"  CI strictly excludes zero? {excludes_zero}")
    log(f"  Verdict: {results[pair]['scientific_verdict']}")
    log("-" * 85)

with open(OUT_JSON, "w") as f:
    json.dump(results, f, indent=2)

with open(OUT_LOG, "w") as f:
    f.write("\n".join(log_lines))

log(f"\nSaved results to: {OUT_JSON}")
log(f"Saved log to: {OUT_LOG}")
