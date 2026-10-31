"""
Compute Genuine Multi-Herpesvirus PSSM (Position-Specific Scoring Matrix) Baseline
==================================================================================
Protocol: AUDIT_PROTOCOL.md Strict Version 4 (Frozen)
Dataset: data/processed/primary_position_folds.tsv (SHA-256 verified)

Builds a true Position-Specific Scoring Matrix (PSSM) from the multi-herpesvirus
multiple sequence alignment:
  - UL23: HSV-1, HSV-2, VZV
  - UL30: HSV-1, HSV-2, VZV, HCMV, HHV-6B

Computes:
  1. Position-specific amino acid probability distribution with Dirichlet / pseudocount smoothing:
     p_pos(a) = (c_pos(a) + alpha * q(a)) / (N_pos + alpha)
     where q(a) is standard background frequency (Robinson & Robinson, 1991).
  2. PSSM Log-Odds Delteriousness Score:
     Score(w -> m) = log2(p_pos(w) / p_pos(m)) + log2(q(m) / q(w))
  3. Contiguous-Block CV (5 blocks) and Position-Grouped CV (5 folds).
  4. Block-level Student's t 95% CI and 99% CI (for contiguous-block means).
  5. Position-cluster bootstrap 95% CI (1,000 resamples over unique positions) for pooled AUROC.
  6. Paired per-block difference against ESM-2 650M (ESM - PSSM).
  7. Paired cluster-bootstrap on pooled AUROC against ESM-2 650M.
"""

import os
import sys
import json
import math
import hashlib
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(ROOT, "scripts"))

from compute_msa_conservation_baseline import orthologs, align_to_ref, assign_block

PRIMARY_TSV = os.path.join(ROOT, "data", "processed", "primary_position_folds.tsv")
SCORED_TSV  = os.path.join(ROOT, "logs", "primary_position_folds_scored.tsv")
PLM_JSON    = os.path.join(ROOT, "logs", "plm_benchmark_results.json")
OUT_JSON    = os.path.join(ROOT, "logs", "pssm_baseline_results.json")
OUT_LOG     = os.path.join(ROOT, "logs", "pssm_baseline.log")

EXPECTED_SHA = "7C4D0A3DE2CD171B07230757F44E0DFBD61C6A0A8570DC8F263195B354C66909"

# Standard background amino acid frequencies (Robinson & Robinson 1991)
BG_FREQ = {
    'A': 0.078, 'R': 0.051, 'N': 0.043, 'D': 0.053, 'C': 0.019,
    'Q': 0.043, 'E': 0.063, 'G': 0.072, 'H': 0.023, 'I': 0.053,
    'L': 0.091, 'K': 0.059, 'M': 0.023, 'F': 0.039, 'P': 0.052,
    'S': 0.068, 'T': 0.059, 'W': 0.014, 'Y': 0.032, 'V': 0.066
}
AMINO_ACIDS = list(BG_FREQ.keys())

# Build reference alignment mapping
pos_residues = {}
for gene, o_dict in orthologs.items():
    ref_seq = o_dict['HSV-1']
    L = len(ref_seq)
    pos_residues[gene] = {p: [ref_seq[p-1]] for p in range(1, L + 1)}
    for o_name, q_seq in o_dict.items():
        if o_name == 'HSV-1':
            continue
        mapping = align_to_ref(ref_seq, q_seq)
        for p in range(1, L + 1):
            ch = mapping.get(p, '-')
            pos_residues[gene][p].append(ch)

def compute_pssm_scores(gene, alpha=1.0):
    """
    Computes PSSM log-odds for all (pos, wt, mut) combinations.
    """
    ref_seq = orthologs[gene]['HSV-1']
    L = len(ref_seq)
    pssm = {}
    
    for p in range(1, L + 1):
        res_list = [ch for ch in pos_residues[gene][p] if ch in BG_FREQ]
        N_p = len(res_list)
        
        # Count observed residues
        counts = {aa: 0 for aa in AMINO_ACIDS}
        for ch in res_list:
            counts[ch] += 1
            
        # Pseudocount-smoothed probabilities
        probs = {}
        for aa in AMINO_ACIDS:
            probs[aa] = (counts[aa] + alpha * BG_FREQ[aa]) / (N_p + alpha)
            
        # Log-odds score relative to background
        log_odds = {aa: math.log2(probs[aa] / BG_FREQ[aa]) for aa in AMINO_ACIDS}
        
        pssm[p] = {
            'probs': probs,
            'log_odds': log_odds
        }
    return pssm

pssm_models = {
    'UL23': compute_pssm_scores('UL23', alpha=1.0),
    'UL30': compute_pssm_scores('UL30', alpha=1.0)
}

def get_pssm_llr(r):
    p = int(r['pos'])
    wt = r['wt']
    mut = r['mut']
    g = r['gene']
    if wt not in BG_FREQ or mut not in BG_FREQ:
        return 0.0
    model = pssm_models[g][p]
    p_wt = model['probs'][wt]
    p_mut = model['probs'][mut]
    # Positive score = WT strongly preferred over mutant by evolution = deleterious/resistant
    return math.log2(p_wt / p_mut)

def get_pssm_log_odds_delta(r):
    p = int(r['pos'])
    wt = r['wt']
    mut = r['mut']
    g = r['gene']
    if wt not in BG_FREQ or mut not in BG_FREQ:
        return 0.0
    model = pssm_models[g][p]
    lo_wt = model['log_odds'][wt]
    lo_mut = model['log_odds'][mut]
    return lo_wt - lo_mut

def main():
    with open(PRIMARY_TSV, "rb") as f:
        act_sha = hashlib.sha256(f.read()).hexdigest().upper()
    if act_sha != EXPECTED_SHA:
        print(f"FATAL: Primary TSV SHA mismatch! {act_sha}")
        sys.exit(1)

    df = pd.read_csv(PRIMARY_TSV, sep='\t')
    df['block_fold_id'] = df.apply(lambda r: assign_block(r['gene'], r['pos']), axis=1)
    df['pssm_llr'] = df.apply(get_pssm_llr, axis=1)

    # Also load scored tsv to obtain ESM-2 650M predictions for pooled bootstrap
    df_scored = pd.read_csv(SCORED_TSV, sep='\t')
    esm_llr_map = dict(zip(zip(df_scored['gene'], df_scored['drug'], df_scored['pos'], df_scored['wt'], df_scored['mut']), df_scored['esm2_650m_llr']))
    df['esm2_650m_score'] = df.apply(lambda r: -esm_llr_map.get((r['gene'], r['drug'], r['pos'], r['wt'], r['mut']), 0.0), axis=1)

    with open(PLM_JSON, "r") as f:
        plm_data = json.load(f)

    pairs = [
        ('UL23', 'Aciclovir', 'UL23-Aciclovir'),
        ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
        ('UL30', 'Foscarnet', 'UL30-Foscarnet')
    ]

    log_lines = []
    def log(msg=""):
        print(msg)
        log_lines.append(msg)

    log("=" * 105)
    log("GENUINE MULTI-HERPESVIRUS PSSM (POSITION-SPECIFIC SCORING MATRIX) EVALUATION")
    log("=" * 105)

    results = {}

    for gene, drug, pair_key in pairs:
        sub = df[(df['gene'] == gene) & (df['drug'] == drug)].copy().reset_index(drop=True)
        y = sub['target_binary'].values.astype(int)
        scores = sub['pssm_llr'].values
        esm_scores = sub['esm2_650m_score'].values
        
        # 1. 5-Fold Position-Grouped CV
        pos_aucs = []
        for f in range(5):
            te_mask = (sub['cv_fold_id'] == f)
            y_te = y[te_mask]
            s_te = scores[te_mask]
            if len(np.unique(y_te)) > 1:
                pos_aucs.append(float(roc_auc_score(y_te, s_te)))
        p_mean = float(np.mean(pos_aucs))
        p_sd   = float(np.std(pos_aucs, ddof=1)) if len(pos_aucs) > 1 else 0.0
        
        # 2. 5-Fold Contiguous-Block CV
        block_aucs = {}
        eval_aucs = []
        for b in range(5):
            te_mask = (sub['block_fold_id'] == b)
            y_te = y[te_mask]
            s_te = scores[te_mask]
            if len(np.unique(y_te)) > 1:
                val = float(roc_auc_score(y_te, s_te))
                block_aucs[f"B{b}"] = val
                eval_aucs.append(val)
            else:
                block_aucs[f"B{b}"] = None
                
        c_mean = float(np.mean(eval_aucs))
        c_sd   = float(np.std(eval_aucs, ddof=1)) if len(eval_aucs) > 1 else 0.0
        n_eval = len(eval_aucs)
        
        # Block Student's t 95% and 99% CI (Strictly for contiguous-block mean)
        df_t = n_eval - 1
        t_crit_95 = float(stats.t.ppf(0.975, df=df_t))
        t_crit_99 = float(stats.t.ppf(0.995, df=df_t))
        se = c_sd / math.sqrt(n_eval)
        t_ci_95 = [max(0.0, c_mean - t_crit_95 * se), min(1.0, c_mean + t_crit_95 * se)]
        t_ci_99 = [max(0.0, c_mean - t_crit_99 * se), min(1.0, c_mean + t_crit_99 * se)]
        
        # 3. Pooled AUROC and Position-cluster bootstrap (Strictly for pooled dataset AUROC)
        pooled_auc = float(roc_auc_score(y, scores))
        pooled_esm_auc = float(roc_auc_score(y, esm_scores))
        
        unique_positions = sub['pos'].unique()
        pos_to_idx = {p: np.where(sub['pos'].values == p)[0] for p in unique_positions}
        np.random.seed(42)
        boot_aucs = []
        boot_deltas = []
        for _ in range(1000):
            boot_pos = np.random.choice(unique_positions, size=len(unique_positions), replace=True)
            boot_idx = np.concatenate([pos_to_idx[p] for p in boot_pos])
            y_boot = y[boot_idx]
            s_boot = scores[boot_idx]
            s_esm_boot = esm_scores[boot_idx]
            if len(np.unique(y_boot)) > 1:
                auc_p = float(roc_auc_score(y_boot, s_boot))
                auc_e = float(roc_auc_score(y_boot, s_esm_boot))
                boot_aucs.append(auc_p)
                boot_deltas.append(auc_e - auc_p)
                
        boot_ci_95 = [float(np.percentile(boot_aucs, 2.5)), float(np.percentile(boot_aucs, 97.5))]
        boot_ci_99 = [float(np.percentile(boot_aucs, 0.5)), float(np.percentile(boot_aucs, 99.5))]
        
        boot_deltas = np.array(boot_deltas)
        paired_boot_p_esm_gt_pssm = float(np.mean(boot_deltas > 0))
        paired_boot_ci_95 = [float(np.percentile(boot_deltas, 2.5)), float(np.percentile(boot_deltas, 97.5))]
        
        # 4. Direct comparison with ESM-2 650M on contiguous blocks
        esm_blocks = plm_data['esm2_650m']['metrics'][pair_key]['per_block_breakdown']
        paired_deltas = []
        esm_vals = []
        pssm_vals = []
        for b_id, val in block_aucs.items():
            if val is not None and esm_blocks[b_id]['auroc'] is not None:
                e_val = esm_blocks[b_id]['auroc']
                paired_deltas.append(e_val - val)
                esm_vals.append(e_val)
                pssm_vals.append(val)
                
        mean_delta = float(np.mean(paired_deltas))
        sd_delta = float(np.std(paired_deltas, ddof=1))
        se_delta = sd_delta / math.sqrt(len(paired_deltas))
        t_stat_d = mean_delta / se_delta if se_delta > 0 else 0.0
        p_val_d = float(2 * (1 - stats.t.cdf(abs(t_stat_d), df=len(paired_deltas) - 1)))
        paired_t_ci_95 = [mean_delta - t_crit_95 * se_delta, mean_delta + t_crit_95 * se_delta]
        
        results[pair_key] = {
            'pos_grouped_mean': p_mean,
            'pos_grouped_sd': p_sd,
            'contig_block_mean': c_mean,
            'contig_block_sd': c_sd,
            'contig_block_evaluable_folds': n_eval,
            'block_student_t_ci_95': t_ci_95,
            'block_student_t_ci_99': t_ci_99,
            'pooled_auroc': pooled_auc,
            'pooled_cluster_bootstrap_ci_95': boot_ci_95,
            'pooled_cluster_bootstrap_ci_99': boot_ci_99,
            'position_cluster_bootstrap_ci_95': boot_ci_95,  # backward compatibility
            'per_block_auroc': block_aucs,
            'esm2_650m_per_block': esm_vals,
            'esm2_650m_mean': float(np.mean(esm_vals)),
            'paired_deltas_esm_minus_pssm': paired_deltas,
            'mean_paired_delta': mean_delta,
            'sd_paired_delta': sd_delta,
            'paired_t_ci_95': paired_t_ci_95,
            'paired_t_stat': t_stat_d,
            'paired_p_val': p_val_d,
            'pooled_esm_auroc': pooled_esm_auc,
            'paired_cluster_boot_delta_mean': float(np.mean(boot_deltas)),
            'paired_cluster_boot_delta_ci_95': paired_boot_ci_95,
            'paired_cluster_boot_p_esm_gt_pssm': paired_boot_p_esm_gt_pssm,
            'passes_stopping_rule_ci95_gt_050': bool(t_ci_95[0] > 0.50),
            'bootstrap_ci95_gt_050': bool(boot_ci_95[0] > 0.50)
        }
        
        log(f"Target Pair: {pair_key}")
        log(f"  Pos-Grouped AUROC:               {p_mean:.4f} +/- {p_sd:.4f}")
        log(f"  Contiguous-Block Mean AUROC:     {c_mean:.4f} +/- {c_sd:.4f}")
        log(f"  Block Student's t 95% CI:        [{t_ci_95[0]:.4f}, {t_ci_95[1]:.4f}]")
        log(f"  Block Student's t 99% CI:        [{t_ci_99[0]:.4f}, {t_ci_99[1]:.4f}]")
        log(f"  Pooled Dataset AUROC:            {pooled_auc:.4f}")
        log(f"  Pooled Cluster Bootstrap 95%:    [{boot_ci_95[0]:.4f}, {boot_ci_95[1]:.4f}]")
        log(f"  Pooled Cluster Bootstrap 99%:    [{boot_ci_99[0]:.4f}, {boot_ci_99[1]:.4f}]")
        log(f"  Per-Block Breakdown:             {block_aucs}")
        log(f"  --- Comparison vs ESM-2 650M ---")
        log(f"  ESM-2 650M Block Mean:           {np.mean(esm_vals):.4f}")
        log(f"  PSSM Block Mean:                 {c_mean:.4f}")
        log(f"  Paired Delta (ESM - PSSM):       {mean_delta:+.4f} +/- {sd_delta:.4f}")
        log(f"  Paired 95% Student's t CI:       [{paired_t_ci_95[0]:.4f}, {paired_t_ci_95[1]:.4f}]")
        log(f"  Paired t-test:                   t = {t_stat_d:.4f}, p = {p_val_d:.4f}")
        log(f"  Pooled Paired Boot Delta Mean:   {np.mean(boot_deltas):+.4f} (95% CI: [{paired_boot_ci_95[0]:.4f}, {paired_boot_ci_95[1]:.4f}])")
        log(f"  Pooled Boot P(ESM > PSSM):       {paired_boot_p_esm_gt_pssm * 100:.1f}%")
        log(f"  Block t-CI > 0.50?               {results[pair_key]['passes_stopping_rule_ci95_gt_050']}")
        log(f"  Bootstrap CI > 0.50?             {results[pair_key]['bootstrap_ci95_gt_050']}")
        log("-" * 105)

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)

    with open(OUT_LOG, "w") as f:
        f.write("\n".join(log_lines))

    log(f"\nSaved PSSM baseline results to: {OUT_JSON}")
    log(f"Saved log to:                  {OUT_LOG}")

if __name__ == '__main__':
    main()
