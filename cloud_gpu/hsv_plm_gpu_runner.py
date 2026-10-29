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
    raise FileNotFoundError("Could not locate primary_position_folds.tsv in candidate paths.")

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
