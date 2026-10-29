"""
HSV-1 Antiviral Resistance Mutation Prediction: Cloud GPU PLM Benchmark Runner
Models:
  1. ESM-2 650M (facebook/esm2_t33_650M_UR50D)
  2. ProtBERT (Rostlab/prot_bert)
  3. ESM-2 35M (facebook/esm2_t12_35M_UR50D) [Baseline]

Protocol:
  - Masked marginal log-likelihood ratio: Delta LLR = log P(mut) - log P(wt)
  - UL23: Single-pass full sequence (376 aa)
  - UL30: Centered-window extraction (max 510 aa + special tokens)
  - Dual-Split Benchmark: Position-Grouped CV (weak) vs Contiguous-Block CV (strict)
  - Stopping Rule: AUROC >= 0.65 with CI > 0.50 under Contiguous-Block on >= 2/3 pairs.
"""

import os
import sys
import json
import time
import math
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from transformers import AutoTokenizer, AutoModelForMaskedLM, BertTokenizer, BertForMaskedLM

REF_SEQS = {
    'UL23': "MASYPCHQHASAFDQAARSRGHNNRRTALRPRRQQKATEVRLEQKMPTLLRVYIDGPHGMGKTTTTQLLVALGSRDDIVYVPEPMTYWRVLGASETIANIYTTQHRLDQGEISAGDAAVVMTSAQITMGMPYAVTDAVLAPHIGGEAGSSHAPPPALTLIFDRHPIAALLCYPAARYLMGSMTPQAVLAFVALIPPTLPGTNIVLGALPEDRHIDRLAKRQRPGERLDLAMLAAIRRVYGLLANTVRYLQGGGSWREDWGQLSGAAVPPQGAEPQSNAGPRPHIGDTLFTLFRAPELLAPNGDLYNVFAWALDVLAKRLRPMHVFILDYDQSPAGCRDALLQLTSGMVQTHVTTPGSIPTICDLARTFAREMGEAN",
    'UL30': "MFSGGGGPLSPGGKSAARAASGFFAPAGPRGASRGPPPCLRQNFYNPYLAPVGTQQKPTGPTQRHTYYSECDEFRFIAPRVLDEDAPPEKRAGVHDGHLKRAPKVYCGGDERDVLRVGSGGFWPRRSRLWGGVDHAPAGFNPTVTVFHVYDILENVEHAYGMRAAQFHARFMDAITPTGTVITLLGLTPEGHRVAVHVYGTRQYFYMNKEEVDRHLQCRAPRDLCERMAAALRESPGASFRGISADHFEAEVVERTDVYYYETRPALFYRVYVRSGRVLSYLCDNFCPAIKKYEGGVDATTRFILDNPGFVTFGWYRLKPGRNNTLAQPRAPMAFGTSSDVEFNCTADNLAIEGGMSDLPAYKLMCFDIECKAGGEDELAFPVAGHPEDLVIQISCLLYDLSTTALEHVLLFSLGSCDLPESHLNELAARGLPTPVVLEFDSEFEMLLAFMTLVKQYGPEFVTGYNIINFDWPFLLAKLTDIYKVPLDGYGRMNGRGVFRVWDIGQSHFQKRSKIKVNGMVNIDMYGIITDKIKLSSYKLNAVAEAVLKDKKKDLSYRDIPAYYAAGPAQRGVIGEYCIQDSLLVGQLFFKFLPHLELSAVARLAGINITRTIYDGQQIRVFTCLLRLADQKGFILPDTQGRFRGAGGEAPKRPAAAREDEERPEEEGEDEDEREEGGGEREPEGARETAGRHVGYQGARVLDPTSGFHVNPVVVFDFASLYPSIIQAHNLCFSTLSLRADAVAHLEAGKDYLEIEVGGRRLFFVKAHVRESLLSILLRDWLAMRKQIRSRIPQSSPEEAVLLDKQQAAIKVVCNSVYGFTGVQHGLLPCLHVAATVTTIGREMLLATREYVHARWAAFEQLLADFPEAADMRAPGPYSMRIIYGDTDSIFVLCRGLTAAGLTAVGDKMASHISRALFLPPIKLECEKTFTKLLLIAKKKYIGVIYGGKMLIKGVDLVRKNNCAFINRTSRALVDLLFYDDTVSGAAAALAERPAEEWLARPLPEGLQAFGAVLVDAHRRITDPERDIQDFVLTAELSRHPRAYTNKRLAHLTVYYKLMARRAQVPSIKDRIPYVIVAQTREVEETVARLAALRELDAAAPGDEPAPPAALPSPAKRPRETPSPADPPGGASKPRKLLVSELAEDPAYAIAHGVALNTDYYFSHLLGAACVTFKALFGNNAKITESLLKRFIPEVWHPPDDVAARLRTAGFGAVGAGATAEETRRMLHRAFDTLA"
}

def assign_contiguous_block(gene, pos):
    if gene == 'UL23':
        block = min(4, (pos - 1) // 75)
    else:
        block = min(4, (pos - 1) // 247)
    return int(block)

EMBEDDED_DATA_GZIP_B64 = "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"

def find_dataset_file():
    candidates = [
        'primary_position_folds.tsv',
        'cloud_gpu/primary_position_folds.tsv',
        'data/processed/primary_position_folds.tsv',
        '../data/processed/primary_position_folds.tsv',
        '/content/primary_position_folds.tsv',
        '/kaggle/input/hsv-primary-folds/primary_position_folds.tsv'
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    
    # Auto-extract embedded dataset if not found on disk
    import gzip, base64
    out_path = 'primary_position_folds.tsv'
    with open(out_path, 'wb') as f:
        f.write(gzip.decompress(base64.b64decode(EMBEDDED_DATA_GZIP_B64)))
    print(f'[Dataset] Extracted self-contained dataset to: {out_path}')
    return out_path

def score_plm_model(model_name, model_type, df, device):
    print(f"\n==================================================================")
    print(f"Loading {model_name} on {device}...")
    print(f"==================================================================")
    
    if model_type == 'esm':
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModelForMaskedLM.from_pretrained(model_name, torch_dtype=torch.float16 if 'cuda' in device else torch.float32)
    elif model_type == 'protbert':
        tokenizer = BertTokenizer.from_pretrained(model_name, do_lower_case=False)
        model = BertForMaskedLM.from_pretrained(model_name, torch_dtype=torch.float16 if 'cuda' in device else torch.float32)
    else:
        raise ValueError(f"Unknown model_type: {model_type}")

    model.to(device)
    model.eval()
    
    canonical_aas = 'ACDEFGHIKLMNPQRSTVWY'
    aa_tokens = {a: tokenizer.convert_tokens_to_ids(a) for a in canonical_aas}
    mask_token_id = tokenizer.mask_token_id

    unique_pos = df[['gene', 'pos']].drop_duplicates().sort_values(['gene', 'pos']).reset_index(drop=True)
    print(f"Scoring {len(unique_pos)} unique sequence positions across UL23 and UL30...")

    pos_llr_cache = {}
    start_time = time.time()

    with torch.no_grad():
        for idx, row in unique_pos.iterrows():
            g = row['gene']
            p = int(row['pos'])
            wt_full = REF_SEQS[g]
            L = len(wt_full)
            p_idx = p - 1  # 0-indexed

            if g == 'UL23':
                # Full sequence single-pass
                window_wt = wt_full
                mask_idx = p_idx
            else:
                # Centered window for UL30 (max 510 aa)
                p_start = max(0, p_idx - 255)
                p_end = min(L, p_start + 510)
                if p_end - p_start < 510 and p_start > 0:
                    p_start = max(0, p_end - 510)
                window_wt = wt_full[p_start:p_end]
                mask_idx = p_idx - p_start

            # Format input
            if model_type == 'protbert':
                # ProtBERT expects space-separated uppercase residues
                tok_list = list(window_wt)
                tok_list[mask_idx] = '[MASK]'
                input_str = ' '.join(tok_list)
                inputs = tokenizer(input_str, return_tensors='pt').to(device)
                seq_mask_pos = (inputs.input_ids[0] == mask_token_id).nonzero(as_tuple=True)[0].item()
            else:
                # ESM-2 expects unspaced sequence with <mask>
                masked_seq = window_wt[:mask_idx] + '<mask>' + window_wt[mask_idx+1:]
                inputs = tokenizer(masked_seq, return_tensors='pt').to(device)
                seq_mask_pos = (inputs.input_ids[0] == mask_token_id).nonzero(as_tuple=True)[0].item()

            logits = model(**inputs).logits[0, seq_mask_pos]
            log_probs = torch.log_softmax(logits, dim=-1)

            wt_aa = wt_full[p_idx]
            wt_tok = aa_tokens.get(wt_aa)
            wt_logp = log_probs[wt_tok].item() if wt_tok is not None else 0.0

            scores = {}
            for aa in canonical_aas:
                tok_id = aa_tokens[aa]
                scores[aa] = float((log_probs[tok_id] - wt_logp).item())

            pos_llr_cache[(g, p)] = scores

            if (idx + 1) % 50 == 0 or (idx + 1) == len(unique_pos):
                elapsed = time.time() - start_time
                print(f"  Processed {idx + 1}/{len(unique_pos)} positions ({elapsed:.1f}s)...", flush=True)

    # Compute Delta LLR for every row in dataset
    llrs = []
    for _, r in df.iterrows():
        g, p, mut = r['gene'], int(r['pos']), r['mut']
        scores = pos_llr_cache.get((g, p), {})
        llrs.append(scores.get(mut, 0.0))
    
    del model
    if 'cuda' in device:
        torch.cuda.empty_cache()

    return llrs

def evaluate_splits(df, score_col):
    pairs = [
        ('UL23', 'Aciclovir', 'UL23-Aciclovir'),
        ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
        ('UL30', 'Foscarnet', 'UL30-Foscarnet')
    ]
    results = {}

    for g, d, p_key in pairs:
        sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
        y = sub['target_binary'].values.astype(int)
        scores = sub[score_col].values

        # Resistance = catalytic disruption (negative Delta LLR), so we invert scores for AUROC
        auc_directed = roc_auc_score(y, -scores)

        # 1. Position-Grouped CV (5 folds)
        pos_fold_aucs = []
        for f in range(5):
            te_mask = (sub['cv_fold_id'] == f)
            y_te = y[te_mask]
            s_te = -scores[te_mask]
            if len(np.unique(y_te)) > 1:
                pos_fold_aucs.append(roc_auc_score(y_te, s_te))
        pos_mean = np.mean(pos_fold_aucs) if pos_fold_aucs else np.nan
        pos_sd = np.std(pos_fold_aucs, ddof=1) if len(pos_fold_aucs) > 1 else 0.0

        # 2. Contiguous-Block CV (5 blocks)
        block_fold_aucs = []
        for b in range(5):
            te_mask = (sub['block_fold_id'] == b)
            y_te = y[te_mask]
            s_te = -scores[te_mask]
            if len(np.unique(y_te)) > 1:
                block_fold_aucs.append(roc_auc_score(y_te, s_te))
        block_mean = np.mean(block_fold_aucs) if block_fold_aucs else np.nan
        block_sd = np.std(block_fold_aucs, ddof=1) if len(block_fold_aucs) > 1 else 0.0

        # 3. Position-Cluster Bootstrap 95% CI under Contiguous-Block holdout
        unique_positions = sub['pos'].unique()
        n_pos = len(unique_positions)
        rng = np.random.RandomState(42)
        b_aucs = []
        for _ in range(500):
            sample_pos = rng.choice(unique_positions, size=n_pos, replace=True)
            boot_idx = [i for p in sample_pos for i in sub[sub['pos'] == p].index]
            y_b = y[boot_idx]
            s_b = -scores[boot_idx]
            if len(np.unique(y_b)) > 1:
                b_aucs.append(roc_auc_score(y_b, s_b))
        
        ci_lower = np.percentile(b_aucs, 2.5) if b_aucs else np.nan
        ci_upper = np.percentile(b_aucs, 97.5) if b_aucs else np.nan

        results[p_key] = {
            'overall_auroc': float(auc_directed),
            'pos_grouped_mean': float(pos_mean),
            'pos_grouped_sd': float(pos_sd),
            'contig_block_mean': float(block_mean),
            'contig_block_sd': float(block_sd),
            'contig_block_ci_95': [float(ci_lower), float(ci_upper)],
            'strictly_passes_stopping_rule': bool(block_mean >= 0.65 and ci_lower > 0.50)
        }

    return results

def main():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"=== HSV-1 ANTIVIRAL RESISTANCE PLM BENCHMARK ===")
    print(f"Device: {device} ({torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'})")

    ds_path = find_dataset_file()
    print(f"Dataset path: {ds_path}")
    df = pd.read_csv(ds_path, sep='\t')
    df['block_fold_id'] = df.apply(lambda r: assign_contiguous_block(r['gene'], r['pos']), axis=1)
    print(f"Loaded {len(df)} primary cohort rows across {df['gene'].nunique()} genes, {df['drug'].nunique()} drugs.")

    models_to_run = [
        ('esm2_35m', 'facebook/esm2_t12_35M_UR50D', 'esm'),
        ('esm2_650m', 'facebook/esm2_t33_650M_UR50D', 'esm'),
        ('protbert', 'Rostlab/prot_bert', 'protbert')
    ]

    all_benchmark_results = {}

    for model_key, model_id, model_type in models_to_run:
        score_col = f'{model_key}_llr'
        try:
            df[score_col] = score_plm_model(model_id, model_type, df, device)
            res = evaluate_splits(df, score_col)
            all_benchmark_results[model_key] = {
                'model_id': model_id,
                'model_type': model_type,
                'metrics': res
            }

            print(f"\n--- {model_id} Results ---")
            print(f"{'Target Pair':<18} | {'Pos-Grouped AUROC':<20} | {'Contig-Block AUROC':<20} | {'95% Bootstrap CI':<20} | {'Passes (>=0.65)?'}")
            print("-" * 95)
            n_passed = 0
            for p_key, m in res.items():
                pass_str = "PASS" if m['strictly_passes_stopping_rule'] else "NO"
                if m['strictly_passes_stopping_rule']:
                    n_passed += 1
                ci_str = f"[{m['contig_block_ci_95'][0]:.4f}, {m['contig_block_ci_95'][1]:.4f}]"
                print(f"{p_key:<18} | {m['pos_grouped_mean']:.4f} +/- {m['pos_grouped_sd']:.4f}   | {m['contig_block_mean']:.4f} +/- {m['contig_block_sd']:.4f}   | {ci_str:<20} | {pass_str}")
            print(f"Stopping rule status: {n_passed}/3 primary pairs passed (Threshold: >= 2/3).")

        except Exception as e:
            print(f"[ERROR] Failed running {model_id}: {str(e)}", file=sys.stderr)
            all_benchmark_results[model_key] = {'error': str(e)}

    # Save output json
    out_json = 'plm_benchmark_results.json'
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(all_benchmark_results, f, indent=2)
    print(f"\nAll benchmark results saved to: {out_json}")

if __name__ == '__main__':
    main()

