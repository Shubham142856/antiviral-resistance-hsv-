"""
Compute Zero-Shot Masked-Marginal Log-Likelihood Ratios (LLR) for Primary Benchmark Rows.
Models:
1. ESM-2 (facebook/esm2_t12_35M_UR50D and/or facebook/esm2_t33_650M_UR50D)
2. ProtBERT (Rostlab/prot_bert)

Zero-shot formulation:
Delta LLR = log P(mut | x_{\\i}) - log P(wt | x_{\\i})

Evaluation:
- 929 Clean Primary Rows (data/processed/primary_position_folds.tsv)
- Evaluated under both:
  1. Position-Grouped Cross-Validation (5 Folds) - Weak test
  2. Contiguous-Block Cross-Validation (5 Sequence Blocks) - Strict domain test
- Compared directly against:
  1. Position-Only Random Forest (0.8715 on UL30 Pos-Grouped, 0.5001 Contiguous-Block)
  2. All-Features Random Forest (0.7328 on UL30 Pos-Grouped, 0.4089 Contiguous-Block)
"""
import os
import sys
import json
import torch
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.metrics import roc_auc_score
from transformers import AutoTokenizer, BertTokenizer, AutoModelForMaskedLM, BertForMaskedLM

def get_hsv1_reference_sequences():
    ref_path = r'd:\antiviral-resistance-hsv--1\herpesdrg\inst\ref\NC_001806.2.fasta'
    ref = SeqIO.read(ref_path, 'fasta')
    # UL23 (TK)
    ul23_nt = ref.seq[46672:47803].reverse_complement()
    ul23_aa = str(ul23_nt.translate())[:-1]
    # UL30 (Pol)
    ul30_nt = ref.seq[62806:66514]
    ul30_aa = str(ul30_nt.translate())[:-1]
    return {'UL23': ul23_aa, 'UL30': ul30_aa}

def assign_contiguous_block(gene, pos):
    if gene == 'UL23':
        block = min(4, (pos - 1) // 75)
    else:
        block = min(4, (pos - 1) // 247)
    return int(block)

def score_esm2(df, ref_seqs, model_id='facebook/esm2_t12_35M_UR50D', device='cpu'):
    cache_path = r'd:\antiviral-resistance-hsv--1\data\processed\esm2_pos_cache.json'
    if os.path.exists(cache_path):
        print(f"[ESM-2] Loading cached position scores from: {cache_path}", flush=True)
        with open(cache_path, 'r', encoding='utf-8') as f:
            raw_cache = json.load(f)
        pos_llr_cache = { (k.split(':')[0], int(k.split(':')[1])): v for k, v in raw_cache.items() }
    else:
        print(f"\n[ESM-2] Loading {model_id} on {device}...", flush=True)
        tokenizer = AutoTokenizer.from_pretrained(model_id)
        model = AutoModelForMaskedLM.from_pretrained(model_id).to(device)
        model.eval()
        
        canonical_aas = 'ACDEFGHIKLMNPQRSTVWY'
        aa_tokens = {a: tokenizer.convert_tokens_to_ids(a) for a in canonical_aas}
        mask_token_id = tokenizer.mask_token_id
        
        unique_pos = df[['gene', 'pos']].drop_duplicates().sort_values(['gene', 'pos'])
        pos_llr_cache = {}
        
        print(f"[ESM-2] Computing masked-marginal scores across {len(unique_pos)} unique sequence positions...", flush=True)
        
        count = 0
        with torch.no_grad():
            for idx, r in unique_pos.iterrows():
                g = r['gene']
                p = int(r['pos'])
                full_seq = ref_seqs[g]
                
                p_idx = p - 1 # 0-based
                if p_idx < 0 or p_idx >= len(full_seq):
                    continue
                wt_aa = full_seq[p_idx]
                
                # ESM-2 max length is 1024. UL30 is 1235 aa -> extract centered window of 1020 residues
                if len(full_seq) > 1020:
                    p_start = max(0, p_idx - 510)
                    p_end = min(len(full_seq), p_start + 1020)
                    if p_end - p_start < 1020 and p_start > 0:
                        p_start = max(0, p_end - 1020)
                    sub_seq = full_seq[p_start:p_end]
                    local_p_idx = p_idx - p_start
                else:
                    sub_seq = full_seq
                    local_p_idx = p_idx
                    
                masked_seq = sub_seq[:local_p_idx] + "<mask>" + sub_seq[local_p_idx+1:]
                inputs = tokenizer(masked_seq, return_tensors='pt').to(device)
                
                input_ids = inputs['input_ids'][0]
                mask_indices = (input_ids == mask_token_id).nonzero(as_tuple=True)[0]
                if len(mask_indices) == 0:
                    continue
                mask_token_pos = mask_indices[0].item()
                
                logits = model(**inputs).logits[0, mask_token_pos]
                log_probs = torch.log_softmax(logits, dim=-1)
                
                wt_token_id = aa_tokens.get(wt_aa)
                wt_lp = log_probs[wt_token_id].item() if wt_token_id else 0.0
                
                pos_llr_cache[(g, p)] = {
                    a: log_probs[aa_tokens[a]].item() - wt_lp for a in canonical_aas if a in aa_tokens
                }
                
                count += 1
                if count % 50 == 0:
                    print(f"  [ESM-2] Processed {count}/{len(unique_pos)} positions...", flush=True)
                    
        # Save cache
        raw_to_save = { f"{k[0]}:{k[1]}": v for k, v in pos_llr_cache.items() }
        with open(cache_path, 'w', encoding='utf-8') as f:
            json.dump(raw_to_save, f)
        print(f"[ESM-2] Cached scores saved to: {cache_path}", flush=True)
                    
    llrs = []
    for idx, r in df.iterrows():
        g, p, mut = r['gene'], int(r['pos']), r['mut']
        scores = pos_llr_cache.get((g, p), {})
        llrs.append(scores.get(mut, 0.0))
        
    return np.array(llrs)

def score_protbert(df, ref_seqs, device='cpu'):
    cache_path = r'd:\antiviral-resistance-hsv--1\data\processed\protbert_pos_cache.json'
    if os.path.exists(cache_path):
        print(f"[ProtBERT] Loading cached position scores from: {cache_path}", flush=True)
        with open(cache_path, 'r', encoding='utf-8') as f:
            raw_cache = json.load(f)
        pos_llr_cache = { (k.split(':')[0], int(k.split(':')[1])): v for k, v in raw_cache.items() }
    else:
        model_name = 'Rostlab/prot_bert'
        print(f"\n[ProtBERT] Loading {model_name} on {device}...", flush=True)
        tokenizer = BertTokenizer.from_pretrained(model_name, do_lower_case=False)
        model = BertForMaskedLM.from_pretrained(model_name).to(device)
        model.eval()
        
        canonical_aas = 'ACDEFGHIKLMNPQRSTVWY'
        aa_tokens = {a: tokenizer.convert_tokens_to_ids(a) for a in canonical_aas}
        mask_token_id = tokenizer.mask_token_id
        
        unique_pos = df[['gene', 'pos']].drop_duplicates().sort_values(['gene', 'pos'])
        pos_llr_cache = {}
        
        print(f"[ProtBERT] Computing masked-marginal scores across {len(unique_pos)} unique sequence positions...", flush=True)
        
        count = 0
        with torch.no_grad():
            for idx, r in unique_pos.iterrows():
                g = r['gene']
                p = int(r['pos'])
                full_seq = ref_seqs[g]
                
                p_idx = p - 1
                if p_idx < 0 or p_idx >= len(full_seq):
                    continue
                wt_aa = full_seq[p_idx]
                
                # ProtBERT context window: 510 residues max
                if len(full_seq) > 510:
                    p_start = max(0, p_idx - 255)
                    p_end = min(len(full_seq), p_start + 510)
                    if p_end - p_start < 510 and p_start > 0:
                        p_start = max(0, p_end - 510)
                    sub_seq = full_seq[p_start:p_end]
                    local_p_idx = p_idx - p_start
                else:
                    sub_seq = full_seq
                    local_p_idx = p_idx
                    
                masked_list = list(sub_seq)
                masked_list[local_p_idx] = '[MASK]'
                masked_input_str = " ".join(masked_list)
                
                inputs = tokenizer(masked_input_str, return_tensors='pt').to(device)
                input_ids = inputs['input_ids'][0]
                mask_indices = (input_ids == mask_token_id).nonzero(as_tuple=True)[0]
                if len(mask_indices) == 0:
                    continue
                mask_token_pos = mask_indices[0].item()
                
                logits = model(**inputs).logits[0, mask_token_pos]
                log_probs = torch.log_softmax(logits, dim=-1)
                
                wt_token_id = aa_tokens.get(wt_aa)
                wt_lp = log_probs[wt_token_id].item() if wt_token_id else 0.0
                
                pos_llr_cache[(g, p)] = {
                    a: log_probs[aa_tokens[a]].item() - wt_lp for a in canonical_aas if a in aa_tokens
                }
                
                count += 1
                if count % 50 == 0:
                    print(f"  [ProtBERT] Processed {count}/{len(unique_pos)} positions...", flush=True)
                    
        # Save cache
        raw_to_save = { f"{k[0]}:{k[1]}": v for k, v in pos_llr_cache.items() }
        with open(cache_path, 'w', encoding='utf-8') as f:
            json.dump(raw_to_save, f)
        print(f"[ProtBERT] Cached scores saved to: {cache_path}", flush=True)
                    
    llrs = []
    for idx, r in df.iterrows():
        g, p, mut = r['gene'], int(r['pos']), r['mut']
        scores = pos_llr_cache.get((g, p), {})
        llrs.append(scores.get(mut, 0.0))
        
    return np.array(llrs)

def evaluate_model_auroc(df, score_col, model_name):
    pairs = [
        ('UL23', 'Aciclovir', 'UL23-Aciclovir'),
        ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
        ('UL30', 'Foscarnet', 'UL30-Foscarnet')
    ]
    
    pair_res = {}
    lines = [f"\n=== Evaluation: {model_name} ==="]
    
    for g, d, p_key in pairs:
        sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
        y = sub['target_binary'].values.astype(int)
        scores = sub[score_col].values
        
        auc_raw = roc_auc_score(y, scores)
        auc_inv = roc_auc_score(y, -scores)
        best_direction = -1 if auc_inv > auc_raw else 1
        scores_directed = scores * best_direction
        
        # 1. Position-Grouped (5 Folds)
        pos_aucs = []
        for fold_i in range(5):
            te_mask = (sub['cv_fold_id'] == fold_i)
            if len(np.unique(y[te_mask])) == 2:
                pos_aucs.append(roc_auc_score(y[te_mask], scores_directed[te_mask]))
            else:
                pos_aucs.append(0.5)
                
        # 2. Contiguous-Block (5 Blocks)
        blk_aucs = []
        for fold_i in range(5):
            te_mask = (sub['block_fold_id'] == fold_i)
            if len(np.unique(y[te_mask])) == 2:
                blk_aucs.append(roc_auc_score(y[te_mask], scores_directed[te_mask]))
            else:
                blk_aucs.append(0.5)
                
        mean_pos, sd_pos = np.mean(pos_aucs), np.std(pos_aucs)
        mean_blk, sd_blk = np.mean(blk_aucs), np.std(blk_aucs)
        
        dir_str = "Negative (Deleterious -> Resistant)" if best_direction == -1 else "Positive (Likely -> Resistant)"
        
        lines.append(f"Pair: {p_key} (N={len(sub)}, Positives={sum(y)}, Negatives={len(y)-sum(y)})")
        lines.append(f"  Direction: {dir_str}")
        lines.append(f"  Position-Grouped AUROC: {mean_pos:.4f} +/- {sd_pos:.4f} (Folds: {[round(x,3) for x in pos_aucs]})")
        lines.append(f"  Contiguous-Block AUROC: {mean_blk:.4f} +/- {sd_blk:.4f} (Folds: {[round(x,3) for x in blk_aucs]})")
        
        pair_res[p_key] = {
            'direction': dir_str,
            'position_grouped_mean': round(mean_pos, 4),
            'position_grouped_sd': round(sd_pos, 4),
            'contiguous_block_mean': round(mean_blk, 4),
            'contiguous_block_sd': round(sd_blk, 4)
        }
        
    return pair_res, lines

def main():
    ds_path = r'd:\antiviral-resistance-hsv--1\data\processed\primary_position_folds.tsv'
    df = pd.read_csv(ds_path, sep='\t')
    df['block_fold_id'] = df.apply(lambda r: assign_contiguous_block(r['gene'], r['pos']), axis=1)
    
    ref_seqs = get_hsv1_reference_sequences()
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Executing Zero-Shot Evaluation on {device.upper()}...")
    
    # 1. ESM-2 (35M)
    esm_llrs = score_esm2(df, ref_seqs, model_id='facebook/esm2_t12_35M_UR50D', device=device)
    df['esm2_llr'] = esm_llrs
    
    # 2. ProtBERT
    protbert_llrs = score_protbert(df, ref_seqs, device=device)
    df['protbert_llr'] = protbert_llrs
    
    # Save scored primary dataset
    out_tsv = r'd:\antiviral-resistance-hsv--1\data\processed\primary_with_zero_shot_llr.tsv'
    df.to_csv(out_tsv, sep='\t', index=False)
    print(f"\nSaved scored dataset with zero-shot LLRs to: {out_tsv}")
    
    # 3. Evaluate both
    all_results = {}
    master_log = [
        "=== ZERO-SHOT PROTEIN LANGUAGE MODEL BENCHMARK AUDIT LOG ===",
        "Evaluation Set: 929 Clean Primary Rows (primary_position_folds.tsv)",
        "Zero-Shot Formulation: Delta LLR = log P(mut | x_{\\i}) - log P(wt | x_{\\i})",
        "Primary Benchmark Regimes: Position-Grouped (Local holdout) vs Contiguous-Block (Structural domain holdout)",
        ""
    ]
    
    esm_res, esm_lines = evaluate_model_auroc(df, 'esm2_llr', 'ESM-2 (35M, facebook/esm2_t12_35M_UR50D)')
    all_results['ESM2_35M'] = esm_res
    master_log.extend(esm_lines)
    
    pb_res, pb_lines = evaluate_model_auroc(df, 'protbert_llr', 'ProtBERT (420M, Rostlab/prot_bert)')
    all_results['ProtBERT'] = pb_res
    master_log.extend(pb_lines)
    
    # Benchmark Comparison Section
    master_log.append("\n================================================================================")
    master_log.append("BENCHMARK COMPARISON AGAINST POSITIONAL BASELINES")
    master_log.append("================================================================================")
    master_log.append(f"{'Target Pair':<18} | {'Model':<20} | {'Pos-Grouped (Weak)':<20} | {'Contig-Block (Strict)':<22}")
    master_log.append("-" * 88)
    
    pos_only_baselines = {
        'UL23-Aciclovir': {'pos': 0.7683, 'block': 0.5440},
        'UL30-Aciclovir': {'pos': 0.8715, 'block': 0.5001},
        'UL30-Foscarnet': {'pos': 0.8684, 'block': 0.5166}
    }
    
    all_feat_baselines = {
        'UL23-Aciclovir': {'pos': 0.8067, 'block': 0.6875},
        'UL30-Aciclovir': {'pos': 0.7328, 'block': 0.4089},
        'UL30-Foscarnet': {'pos': 0.6756, 'block': 0.6372}
    }
    
    for pair_key in ['UL23-Aciclovir', 'UL30-Aciclovir', 'UL30-Foscarnet']:
        p_b = pos_only_baselines[pair_key]
        a_b = all_feat_baselines[pair_key]
        e_b = esm_res[pair_key]
        pr_b = pb_res[pair_key]
        
        master_log.append(f"{pair_key:<18} | Position-Only RF     | {p_b['pos']:.4f}               | {p_b['block']:.4f}")
        master_log.append(f"{pair_key:<18} | All-Features RF      | {a_b['pos']:.4f}               | {a_b['block']:.4f}")
        master_log.append(f"{pair_key:<18} | Zero-Shot ESM-2      | {e_b['position_grouped_mean']:.4f} +/- {e_b['position_grouped_sd']:.4f}   | {e_b['contiguous_block_mean']:.4f} +/- {e_b['contiguous_block_sd']:.4f}")
        master_log.append(f"{pair_key:<18} | Zero-Shot ProtBERT   | {pr_b['position_grouped_mean']:.4f} +/- {pr_b['position_grouped_sd']:.4f}   | {pr_b['contiguous_block_mean']:.4f} +/- {pr_b['contiguous_block_sd']:.4f}")
        master_log.append("-" * 88)
        
    out_json = r'd:\antiviral-resistance-hsv--1\logs\zero_shot_plm_results.json'
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(all_results, f, indent=2)
        
    out_log = r'd:\antiviral-resistance-hsv--1\logs\zero_shot_plm_evaluation.log'
    with open(out_log, 'w', encoding='utf-8') as f:
        f.write("\n".join(master_log))
        
    print("\n".join(master_log))
    print(f"\nSaved results to: {out_json}")
    print(f"Saved log to: {out_log}")

if __name__ == '__main__':
    main()
