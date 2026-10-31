"""
Paired Dissociation Tests: ESM-2 650M vs Alignment Baselines
=============================================================
Computes:
  1. ESM-2 650M vs Negative MSA Positional Entropy (saved to paired_dissociation_entropy_test.json)
  2. ESM-2 650M vs Strongest Alignment Baseline per target pair (saved to paired_dissociation_strongest_test.json):
     - UL23-Aciclovir: Multi-Herpesvirus Weighted Conservation (0.7373)
     - UL30-Aciclovir: Negative MSA Positional Entropy (0.6999)
     - UL30-Foscarnet: Negative MSA Positional Entropy (0.8246)
"""

import os
import json
import math
import numpy as np
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PLM_JSON = os.path.join(ROOT, "logs", "plm_benchmark_results.json")
MSA_JSON = os.path.join(ROOT, "logs", "msa_conservation_results.json")

OUT_ENTROPY_JSON = os.path.join(ROOT, "logs", "paired_dissociation_entropy_test.json")
OUT_ENTROPY_LOG  = os.path.join(ROOT, "logs", "paired_dissociation_entropy_test.log")

OUT_STRONGEST_JSON = os.path.join(ROOT, "logs", "paired_dissociation_strongest_test.json")
OUT_STRONGEST_LOG  = os.path.join(ROOT, "logs", "paired_dissociation_strongest_test.log")

with open(PLM_JSON, "r") as f:
    plm_data = json.load(f)

with open(MSA_JSON, "r") as f:
    msa_data = json.load(f)

pairs = ['UL23-Aciclovir', 'UL30-Aciclovir', 'UL30-Foscarnet']

# -------------------------------------------------------------------------
# Test 1: ESM-2 650M vs Negative Positional Entropy (Original Protocol)
# -------------------------------------------------------------------------
entropy_results = {}
entropy_log_lines = []

def log_entropy(msg=""):
    print(msg)
    entropy_log_lines.append(msg)

log_entropy("=" * 85)
log_entropy("PAIRED DISSOCIATION TEST: ZERO-SHOT ESM-2 650M vs POSITIONAL ENTROPY")
log_entropy("=" * 85)

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
    
    entropy_results[pair] = {
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
    
    log_entropy(f"Target Pair: {pair}")
    log_entropy(f"  Evaluable Blocks: {block_names} (N={n_eval})")
    log_entropy(f"  ESM-2 650M per-block: {[round(x, 4) for x in esm_vals]} (Mean = {np.mean(esm_vals):.4f})")
    log_entropy(f"  Entropy    per-block: {[round(x, 4) for x in msa_vals]} (Mean = {np.mean(msa_vals):.4f})")
    log_entropy(f"  Paired deltas (ESM - Entropy): {[round(x, 4) for x in paired_deltas]}")
    log_entropy(f"  Mean paired delta: {mean_delta:+.4f} +/- {sd_delta:.4f} (SE = {se:.4f})")
    log_entropy(f"  t({df_t}) = {t_stat:.4f}, p = {p_val:.4f}")
    log_entropy(f"  95% Student's t CI: [{ci_95[0]:.4f}, {ci_95[1]:.4f}]")
    log_entropy(f"  CI strictly excludes zero? {excludes_zero}")
    log_entropy(f"  Verdict: {entropy_results[pair]['scientific_verdict']}")
    log_entropy("-" * 85)

with open(OUT_ENTROPY_JSON, "w") as f:
    json.dump(entropy_results, f, indent=2)

with open(OUT_ENTROPY_LOG, "w") as f:
    f.write("\n".join(entropy_log_lines))

# -------------------------------------------------------------------------
# Test 2: ESM-2 650M vs Strongest Alignment Baseline (Table 3 Exact Match)
# -------------------------------------------------------------------------
strongest_spec = {
    'UL23-Aciclovir': 'Multi-Herpesvirus Weighted Conservation',
    'UL30-Aciclovir': 'Negative MSA Entropy (-Entropy)',
    'UL30-Foscarnet': 'Negative MSA Entropy (-Entropy)'
}

strongest_results = {}
strongest_log_lines = []

def log_strongest(msg=""):
    print(msg)
    strongest_log_lines.append(msg)

log_strongest("=" * 85)
log_strongest("PAIRED DISSOCIATION TEST: ZERO-SHOT ESM-2 650M vs STRONGEST BASELINE (TABLE 3)")
log_strongest("=" * 85)

for pair in pairs:
    comparator_name = strongest_spec[pair]
    esm_blocks = plm_data['esm2_650m']['metrics'][pair]['per_block_breakdown']
    comp_blocks = msa_data['results_by_pair'][pair][comparator_name]['per_block_auroc']
    
    paired_deltas = []
    block_names = []
    esm_vals = []
    comp_vals = []
    
    for b_id in ['B0', 'B1', 'B2', 'B3', 'B4']:
        e_val = esm_blocks[b_id]['auroc']
        c_val = comp_blocks.get(b_id)
        
        if e_val is not None and c_val is not None:
            block_names.append(b_id)
            esm_vals.append(e_val)
            comp_vals.append(c_val)
            delta = e_val - c_val
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
    
    strongest_results[pair] = {
        'strongest_baseline_name': comparator_name,
        'evaluable_blocks': block_names,
        'n_eval': n_eval,
        'esm2_650m_per_block': esm_vals,
        'strongest_baseline_per_block': comp_vals,
        'esm2_650m_mean': float(np.mean(esm_vals)),
        'strongest_baseline_mean': float(np.mean(comp_vals)),
        'paired_deltas': paired_deltas,
        'mean_paired_delta': mean_delta,
        'sd_paired_delta': sd_delta,
        'se_paired_delta': se,
        't_statistic': t_stat,
        'degrees_of_freedom': df_t,
        'p_value': p_val,
        'student_t_ci_95': ci_95,
        'strictly_excludes_zero': excludes_zero,
        'table3_label': 'tied' if p_val > 0.50 and abs(mean_delta) < 0.01 else ('n.s.' if not excludes_zero else 'sig')
    }
    
    log_strongest(f"Target Pair: {pair}")
    log_strongest(f"  Strongest Baseline: {comparator_name}")
    log_strongest(f"  Evaluable Blocks: {block_names} (N={n_eval})")
    log_strongest(f"  ESM-2 650M per-block: {[round(x, 4) for x in esm_vals]} (Mean = {np.mean(esm_vals):.4f})")
    log_strongest(f"  Baseline   per-block: {[round(x, 4) for x in comp_vals]} (Mean = {np.mean(comp_vals):.4f})")
    log_strongest(f"  Paired deltas: {[round(x, 4) for x in paired_deltas]}")
    log_strongest(f"  Mean paired delta: {mean_delta:+.4f} +/- {sd_delta:.4f} (SE = {se:.4f})")
    log_strongest(f"  t({df_t}) = {t_stat:.4f}, p = {p_val:.4f}")
    log_strongest(f"  95% Student's t CI: [{ci_95[0]:.4f}, {ci_95[1]:.4f}]")
    log_strongest(f"  Table 3 Designation: {strongest_results[pair]['table3_label']}")
    log_strongest("-" * 85)

with open(OUT_STRONGEST_JSON, "w") as f:
    json.dump(strongest_results, f, indent=2)

with open(OUT_STRONGEST_LOG, "w") as f:
    f.write("\n".join(strongest_log_lines))

print(f"\nSaved Entropy test results to: {OUT_ENTROPY_JSON}")
print(f"Saved Strongest test results to: {OUT_STRONGEST_JSON}")
