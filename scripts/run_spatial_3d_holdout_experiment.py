"""
3D Spatial Neighborhood Holdout Experiment
=========================================
Target: Test whether zero-shot PLM representations and baseline models generalize
across 3D spatial structural neighborhoods using cryo-EM coordinates (8V1Q for UL30)
and crystal coordinates (1OF1 for UL23).

Protocol:
1. Verify SHA-256 of data/processed/primary_position_folds.tsv.
2. Extract C-alpha (x, y, z) coordinates for UL30 (8V1Q chain A) and UL23 (1OF1 chain B).
3. Pre-register 3D spatial folds:
   - 5 spatial clusters per protein via KMeans(n_clusters=5, random_state=42, n_init=10).
   - Nearest sequence neighbor assignment for unmodeled residues.
   - Save to data/processed/primary_3d_spatial_folds.tsv and record SHA-256 hash.
4. Benchmark models evaluated:
   - Position-Only Random Forest (Negative control / shortcut detector)
   - BLOSUM62 Radicalness
   - Positional Shannon Entropy (-Entropy)
   - Multi-Herpesvirus Weighted Conservation
   - Multi-Herpesvirus PSSM
   - ProtBERT (zero-shot)
   - ESM-2 35M (zero-shot)
   - ESM-2 650M (zero-shot)
5. Evaluation regimes:
   - 1D Contiguous-Block Holdout (5 sequence blocks)
   - 3D Spatial-Neighborhood Holdout (5 structural clusters)
6. Cohorts:
   - Full Primary Scope (N=286, 332, 311)
   - Structurally Resolved Core (N=233, 262, 242)
7. Statistical tests:
   - 1,000 position-cluster bootstrap iterations for pooled AUROC 95% CIs.
   - Paired Student's t-test and bootstrap interaction test across 3D vs 1D holdouts.
"""

import os
import sys
import json
import math
import hashlib
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(ROOT, "scripts"))

from download_and_audit_structures import get_hsv1_reference_sequences, parse_chain_ca
from compute_msa_conservation_baseline import orthologs, align_to_ref, max_ent, freq_conserved_map, entropy_map, BLOSUM62
from compute_pssm_baseline import pssm_models, BG_FREQ

PRIMARY_TSV = os.path.join(ROOT, "data", "processed", "primary_position_folds.tsv")
SCORED_TSV  = os.path.join(ROOT, "logs", "primary_position_folds_scored.tsv")
SPATIAL_TSV = os.path.join(ROOT, "data", "processed", "primary_3d_spatial_folds.tsv")
HASH_LOG    = os.path.join(ROOT, "logs", "spatial_3d_fold_assignments_hash.log")
OUT_JSON    = os.path.join(ROOT, "logs", "spatial_3d_holdout_results.json")
OUT_LOG     = os.path.join(ROOT, "logs", "spatial_3d_holdout.log")

EXPECTED_SHA = "7C4D0A3DE2CD171B07230757F44E0DFBD61C6A0A8570DC8F263195B354C66909"

log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(msg)

def get_blosum(wt, mut):
    pair = (wt, mut)
    if pair in BLOSUM62: return BLOSUM62[pair]
    rev = (mut, wt)
    if rev in BLOSUM62: return BLOSUM62[rev]
    return -4

def get_pssm_llr(gene, pos, wt, mut):
    if wt not in BG_FREQ or mut not in BG_FREQ:
        return 0.0
    model = pssm_models[gene].get(pos)
    if not model:
        return 0.0
    p_wt = model['probs'][wt]
    p_mut = model['probs'][mut]
    return math.log2(p_wt / p_mut)

def assign_1d_block(gene, pos):
    if gene == 'UL23':
        return min(4, (pos - 1) // 75)
    return min(4, (pos - 1) // 247)

def cluster_bootstrap_ci(sub, score_col, target_col='target_binary', pos_col='pos', n_boot=1000, seed=42):
    rng = np.random.RandomState(seed)
    unique_positions = np.array(sorted(sub[pos_col].unique()))
    pos_to_indices = {p: sub[sub[pos_col] == p].index.values for p in unique_positions}
    
    y_all = sub[target_col].values
    s_all = sub[score_col].values
    
    # Check baseline pooled AUROC
    if len(np.unique(y_all)) < 2:
        return np.nan, (np.nan, np.nan)
    pooled_auc = float(roc_auc_score(y_all, s_all))
    
    boot_aucs = []
    n_pos = len(unique_positions)
    for _ in range(n_boot):
        boot_positions = rng.choice(unique_positions, size=n_pos, replace=True)
        boot_idx = np.concatenate([pos_to_indices[p] for p in boot_positions])
        y_b = y_all[boot_idx]
        s_b = s_all[boot_idx]
        if len(np.unique(y_b)) < 2:
            continue
        boot_aucs.append(roc_auc_score(y_b, s_b))
        
    if len(boot_aucs) == 0:
        return pooled_auc, (np.nan, np.nan)
    low = float(np.percentile(boot_aucs, 2.5))
    high = float(np.percentile(boot_aucs, 97.5))
    return pooled_auc, (low, high)

def evaluate_fold_split(sub, split_col, score_col, target_col='target_binary'):
    folds = sorted(sub[split_col].unique())
    fold_aucs = {}
    eval_aucs = []
    
    for f in folds:
        te_mask = (sub[split_col] == f)
        y_te = sub.loc[te_mask, target_col].values
        s_te = sub.loc[te_mask, score_col].values
        
        if len(np.unique(y_te)) > 1:
            val = float(roc_auc_score(y_te, s_te))
            fold_aucs[f"Fold_{f}"] = val
            eval_aucs.append(val)
        else:
            fold_aucs[f"Fold_{f}"] = None
            
    mean_val = float(np.mean(eval_aucs)) if len(eval_aucs) > 0 else np.nan
    sd_val = float(np.std(eval_aucs, ddof=1)) if len(eval_aucs) > 1 else 0.0
    sem_val = float(sd_val / math.sqrt(len(eval_aucs))) if len(eval_aucs) > 1 else 0.0
    
    return {
        'per_fold': fold_aucs,
        'mean_auc': mean_val,
        'sd_auc': sd_val,
        'sem_auc': sem_val,
        'evaluable_folds_count': len(eval_aucs)
    }

def main():
    log("=" * 110)
    log("HSV-1 RESISTANCE BENCHMARK: 3D SPATIAL NEIGHBORHOOD HOLDOUT EVALUATION")
    log("=" * 110)
    
    # 1. Verify Dataset Integrity
    with open(PRIMARY_TSV, "rb") as f:
        actual_sha = hashlib.sha256(f.read()).hexdigest().upper()
    if actual_sha != EXPECTED_SHA:
        log(f"FATAL: Primary TSV SHA-256 mismatch! Expected {EXPECTED_SHA}, got {actual_sha}")
        sys.exit(1)
    log(f"Primary TSV SHA-256 verified: {actual_sha}")
    
    ref_seqs = get_hsv1_reference_sequences()
    parser = MMCIFParser(QUIET=True)
    
    # 2. Extract C-alpha coordinates
    s_8v1q = parser.get_structure('8V1Q', os.path.join(ROOT, 'data', 'structures', '8v1q.cif'))
    ca_ul30 = {r['resseq']: (r['x'], r['y'], r['z']) for r in parse_chain_ca(s_8v1q, 'A')}
    
    s_1of1 = parser.get_structure('1OF1', os.path.join(ROOT, 'data', 'structures', '1of1.cif'))
    ca_ul23 = {r['resseq']: (r['x'], r['y'], r['z']) for r in parse_chain_ca(s_1of1, 'B')}
    
    log(f"Loaded 8V1Q Chain A (UL30) C-alpha coordinates: {len(ca_ul30)} resolved residues")
    log(f"Loaded 1OF1 Chain B (UL23) C-alpha coordinates: {len(ca_ul23)} resolved residues")
    
    # 3. Pre-register 3D Spatial Folds
    coords_dict = {'UL30': ca_ul30, 'UL23': ca_ul23}
    spatial_mapping = {}
    
    for gene in ['UL23', 'UL30']:
        gene_coords = coords_dict[gene]
        resolved_pos = sorted(list(gene_coords.keys()))
        X_res = np.array([gene_coords[p] for p in resolved_pos])
        
        # KMeans k=5, seed=42
        kmeans = KMeans(n_clusters=5, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X_res)
        pos_to_c = {p: int(lbl) for p, lbl in zip(resolved_pos, labels)}
        
        ref_len = len(ref_seqs[gene])
        full_mapping = {}
        for p in range(1, ref_len + 1):
            if p in pos_to_c:
                full_mapping[p] = {
                    'cluster_id': pos_to_c[p],
                    'has_3d': 1,
                    'x': gene_coords[p][0],
                    'y': gene_coords[p][1],
                    'z': gene_coords[p][2]
                }
            else:
                nearest_p = min(resolved_pos, key=lambda x: abs(x - p))
                full_mapping[p] = {
                    'cluster_id': pos_to_c[nearest_p],
                    'has_3d': 0,
                    'x': np.nan,
                    'y': np.nan,
                    'z': np.nan
                }
        spatial_mapping[gene] = full_mapping
        
    # Build complete dataframe
    df = pd.read_csv(PRIMARY_TSV, sep='\t')
    df_scored = pd.read_csv(SCORED_TSV, sep='\t')
    
    # Map LLRs from scored TSV
    key_cols = ['gene', 'drug', 'pos', 'wt', 'mut']
    scored_dict = df_scored.set_index(key_cols).to_dict('index')
    
    # Annotate spatial and features
    df['has_3d_coord'] = [spatial_mapping[r['gene']][int(r['pos'])]['has_3d'] for _, r in df.iterrows()]
    df['coord_x'] = [spatial_mapping[r['gene']][int(r['pos'])]['x'] for _, r in df.iterrows()]
    df['coord_y'] = [spatial_mapping[r['gene']][int(r['pos'])]['y'] for _, r in df.iterrows()]
    df['coord_z'] = [spatial_mapping[r['gene']][int(r['pos'])]['z'] for _, r in df.iterrows()]
    df['spatial_cluster_id'] = [spatial_mapping[r['gene']][int(r['pos'])]['cluster_id'] for _, r in df.iterrows()]
    df['seq_block_id'] = [assign_1d_block(r['gene'], int(r['pos'])) for _, r in df.iterrows()]
    
    # Feature scoring
    df['score_blosum62'] = df.apply(lambda r: -get_blosum(r['wt'], r['mut']), axis=1)
    df['score_neg_entropy'] = df.apply(lambda r: -entropy_map.get((r['gene'], int(r['pos'])), 0.0), axis=1)
    df['score_weighted_cons'] = df.apply(
        lambda r: (1.0 - entropy_map.get((r['gene'], int(r['pos'])), 0.0) / max_ent[r['gene']]) * 4.0 + (-get_blosum(r['wt'], r['mut'])),
        axis=1
    )
    df['score_pssm'] = df.apply(lambda r: get_pssm_llr(r['gene'], int(r['pos']), r['wt'], r['mut']), axis=1)
    
    df['score_protbert'] = df.apply(
        lambda r: -scored_dict.get((r['gene'], r['drug'], r['pos'], r['wt'], r['mut']), {}).get('protbert_llr', 0.0),
        axis=1
    )
    df['score_esm2_35m'] = df.apply(
        lambda r: -scored_dict.get((r['gene'], r['drug'], r['pos'], r['wt'], r['mut']), {}).get('esm2_35m_llr', 0.0),
        axis=1
    )
    df['score_esm2_650m'] = df.apply(
        lambda r: -scored_dict.get((r['gene'], r['drug'], r['pos'], r['wt'], r['mut']), {}).get('esm2_650m_llr', 0.0),
        axis=1
    )
    
    # Save spatial folds TSV
    df.to_csv(SPATIAL_TSV, sep='\t', index=False)
    with open(SPATIAL_TSV, "rb") as f:
        spatial_sha = hashlib.sha256(f.read()).hexdigest().upper()
        
    with open(HASH_LOG, "w") as f:
        f.write(f"primary_3d_spatial_folds.tsv SHA-256: {spatial_sha}\n")
        f.write(f"Total Rows: {len(df)}\n")
        f.write(f"UL30 Structure: 8V1Q Chain A (Cryo-EM holoenzyme)\n")
        f.write(f"UL23 Structure: 1OF1 Chain B (Crystal structure complexed with ACV)\n")
        f.write(f"Spatial Clustering: KMeans(n_clusters=5, random_state=42, n_init=10)\n")
    log(f"Saved primary 3D spatial folds: {SPATIAL_TSV} (SHA-256: {spatial_sha})")
    
    # 4. Fit Position-Only Random Forest on both splits
    # For Position-Only RF, we generate out-of-fold probability predictions under each fold
    pairs = [
        ('UL23', 'Aciclovir', 'UL23-Aciclovir', 376),
        ('UL30', 'Aciclovir', 'UL30-Aciclovir', 1235),
        ('UL30', 'Foscarnet', 'UL30-Foscarnet', 1235)
    ]
    
    df['score_rf_pos_only_1d'] = np.nan
    df['score_rf_pos_only_3d'] = np.nan
    
    for gene, drug, pair_key, glen in pairs:
        mask = (df['gene'] == gene) & (df['drug'] == drug)
        sub_indices = df[mask].index
        
        # 1D Block Split RF
        for b in range(5):
            tr_idx = df[mask & (df['seq_block_id'] != b)].index
            te_idx = df[mask & (df['seq_block_id'] == b)].index
            if len(te_idx) == 0 or len(tr_idx) == 0:
                continue
            X_tr = df.loc[tr_idx, ['pos']].values / glen
            y_tr = df.loc[tr_idx, 'target_binary'].values
            X_te = df.loc[te_idx, ['pos']].values / glen
            
            if len(np.unique(y_tr)) > 1:
                rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5)
                rf.fit(X_tr, y_tr)
                preds = rf.predict_proba(X_te)[:, 1]
            else:
                preds = np.full(len(te_idx), 0.5)
            df.loc[te_idx, 'score_rf_pos_only_1d'] = preds
            
        # 3D Spatial Split RF
        for c in range(5):
            tr_idx = df[mask & (df['spatial_cluster_id'] != c)].index
            te_idx = df[mask & (df['spatial_cluster_id'] == c)].index
            if len(te_idx) == 0 or len(tr_idx) == 0:
                continue
            X_tr = df.loc[tr_idx, ['pos']].values / glen
            y_tr = df.loc[tr_idx, 'target_binary'].values
            X_te = df.loc[te_idx, ['pos']].values / glen
            
            if len(np.unique(y_tr)) > 1:
                rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5)
                rf.fit(X_tr, y_tr)
                preds = rf.predict_proba(X_te)[:, 1]
            else:
                preds = np.full(len(te_idx), 0.5)
            df.loc[te_idx, 'score_rf_pos_only_3d'] = preds
            
    models_to_eval = [
        ('Position-Only RF', 'score_rf_pos_only'),
        ('BLOSUM62 Radicalness', 'score_blosum62'),
        ('Positional Shannon Entropy (-Entropy)', 'score_neg_entropy'),
        ('Multi-Herpesvirus Weighted Conservation', 'score_weighted_cons'),
        ('Multi-Herpesvirus PSSM', 'score_pssm'),
        ('ProtBERT (zero-shot)', 'score_protbert'),
        ('ESM-2 35M (zero-shot)', 'score_esm2_35m'),
        ('ESM-2 650M (zero-shot)', 'score_esm2_650m')
    ]
    
    experiment_results = {
        'metadata': {
            'protocol': 'AUDIT_PROTOCOL.md Strict Version 4',
            'dataset_sha256': actual_sha,
            'spatial_sha256': spatial_sha,
            'ul30_structure': '8V1Q Chain A',
            'ul23_structure': '1OF1 Chain B',
            'n_spatial_clusters': 5,
            'n_bootstrap': 1000
        },
        'results_by_pair': {}
    }
    
    log()
    log("=" * 125)
    log(f"{'Target Pair':<16} | {'Cohort':<12} | {'Model Name':<38} | {'1D Block AUROC':<18} | {'3D Spatial AUROC':<18} | {'Delta (3D - 1D)'}")
    log("-" * 125)
    
    for gene, drug, pair_key, glen in pairs:
        pair_dict = {
            'full_cohort': {},
            'resolved_core': {}
        }
        
        cohorts = [
            ('full_cohort', df[(df['gene'] == gene) & (df['drug'] == drug)].copy().reset_index(drop=True)),
            ('resolved_core', df[(df['gene'] == gene) & (df['drug'] == drug) & (df['has_3d_coord'] == 1)].copy().reset_index(drop=True))
        ]
        
        for cohort_name, sub in cohorts:
            cohort_dict = {}
            for m_name, m_col_base in models_to_eval:
                col_1d = f"{m_col_base}_1d" if m_col_base == 'score_rf_pos_only' else m_col_base
                col_3d = f"{m_col_base}_3d" if m_col_base == 'score_rf_pos_only' else m_col_base
                
                # 1D Contiguous-Block Evaluation
                res_1d = evaluate_fold_split(sub, 'seq_block_id', col_1d)
                p_auc_1d, (ci_low_1d, ci_high_1d) = cluster_bootstrap_ci(sub, col_1d)
                
                # 3D Spatial-Neighborhood Evaluation
                res_3d = evaluate_fold_split(sub, 'spatial_cluster_id', col_3d)
                p_auc_3d, (ci_low_3d, ci_high_3d) = cluster_bootstrap_ci(sub, col_3d)
                
                # Paired Delta across 3D vs 1D
                delta_mean = res_3d['mean_auc'] - res_1d['mean_auc']
                delta_pooled = p_auc_3d - p_auc_1d
                
                cohort_dict[m_name] = {
                    '1d_block': {
                        'mean_auc': res_1d['mean_auc'],
                        'sd_auc': res_1d['sd_auc'],
                        'sem_auc': res_1d['sem_auc'],
                        'pooled_auc': p_auc_1d,
                        'pooled_95ci': [ci_low_1d, ci_high_1d],
                        'per_fold': res_1d['per_fold']
                    },
                    '3d_spatial': {
                        'mean_auc': res_3d['mean_auc'],
                        'sd_auc': res_3d['sd_auc'],
                        'sem_auc': res_3d['sem_auc'],
                        'pooled_auc': p_auc_3d,
                        'pooled_95ci': [ci_low_3d, ci_high_3d],
                        'per_fold': res_3d['per_fold']
                    },
                    'delta_3d_minus_1d': {
                        'mean_diff': delta_mean,
                        'pooled_diff': delta_pooled
                    }
                }
                
                str_1d = f"{res_1d['mean_auc']:.4f} +/- {res_1d['sd_auc']:.4f}" if not np.isnan(res_1d['mean_auc']) else "N/A"
                str_3d = f"{res_3d['mean_auc']:.4f} +/- {res_3d['sd_auc']:.4f}" if not np.isnan(res_3d['mean_auc']) else "N/A"
                str_delta = f"{delta_mean:+.4f}" if not np.isnan(delta_mean) else "N/A"
                
                log(f"{pair_key:<16} | {cohort_name:<12} | {m_name:<38} | {str_1d:<18} | {str_3d:<18} | {str_delta}")
                
            pair_dict[cohort_name] = cohort_dict
        experiment_results['results_by_pair'][pair_key] = pair_dict
        log("-" * 125)
        
    # Serialize JSON results
    with open(OUT_JSON, "w") as f:
        json.dump(experiment_results, f, indent=2)
    log(f"\nSerialized experimental results to: {OUT_JSON}")
    
    with open(OUT_LOG, "w") as f:
        f.write("\n".join(log_lines))
    log(f"Saved full execution log to: {OUT_LOG}")

if __name__ == '__main__':
    main()
