import pathlib
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

def assign_contiguous_block(gene, pos):
    if gene == 'UL23':
        block = min(4, (pos - 1) // 75)
    else:
        block = min(4, (pos - 1) // 247)
    return int(block)

def main():
    ds_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/primary_position_folds.tsv')
    cache_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/esm2_pos_cache.json')
    
    df = pd.read_csv(ds_path, sep='\t')
    df['block_fold_id'] = df.apply(lambda r: assign_contiguous_block(r['gene'], r['pos']), axis=1)
    
    with open(cache_path, 'r', encoding='utf-8') as f:
        raw_cache = json.load(f)
    pos_llr_cache = { (k.split(':')[0], int(k.split(':')[1])): v for k, v in raw_cache.items() }
    
    llrs = []
    for idx, r in df.iterrows():
        g, p, mut = r['gene'], int(r['pos']), r['mut']
        scores = pos_llr_cache.get((g, p), {})
        llrs.append(scores.get(mut, 0.0))
    df['esm2_llr'] = llrs
    
    pairs = [
        ('UL23', 'Aciclovir', 'UL23-Aciclovir'),
        ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
        ('UL30', 'Foscarnet', 'UL30-Foscarnet')
    ]
    
    print("\n=== ZERO-SHOT ESM-2 (35M) EVALUATION RESULTS ===")
    print(f"{'Target Pair':<18} | {'Direction':<22} | {'Pos-Grouped AUROC':<22} | {'Contig-Block AUROC':<22}")
    print("-" * 90)
    
    for g, d, p_key in pairs:
        sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
        y = sub['target_binary'].values.astype(int)
        scores = sub['esm2_llr'].values
        
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
        
        dir_str = "Deleterious -> Res" if best_direction == -1 else "Benign -> Res"
        print(f"{p_key:<18} | {dir_str:<22} | {mean_pos:.4f} +/- {sd_pos:.4f}       | {mean_blk:.4f} +/- {sd_blk:.4f}")
        print(f"   Pos folds: {[round(x, 3) for x in pos_aucs]}")
        print(f"   Blk folds: {[round(x, 3) for x in blk_aucs]}")

if __name__ == '__main__':
    main()

