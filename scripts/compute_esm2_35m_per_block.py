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

df = pd.read_csv(str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/primary_position_folds.tsv'), sep='\t')
df['block'] = df.apply(lambda r: assign_contiguous_block(r['gene'], r['pos']), axis=1)

with open(str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/esm2_pos_cache.json'), 'r', encoding='utf-8') as f:
    raw_cache = json.load(f)
pos_llr_cache = { (k.split(':')[0], int(k.split(':')[1])): v for k, v in raw_cache.items() }

df['esm2_35m_llr'] = [pos_llr_cache.get((r['gene'], int(r['pos'])), {}).get(r['mut'], 0.0) for _, r in df.iterrows()]

pairs = [
    ('UL23', 'Aciclovir', 'UL23-Aciclovir'),
    ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
    ('UL30', 'Foscarnet', 'UL30-Foscarnet')
]

print("=== ESM-2 (35M) PER-BLOCK EVALUATION ===")
for g, d, name in pairs:
    sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
    y = sub['target_binary'].values.astype(int)
    scores = -sub['esm2_35m_llr'].values
    
    b_aucs = []
    b_strs = []
    for b in range(5):
        m = sub['block'] == b
        y_b = y[m]
        s_b = scores[m]
        if len(np.unique(y_b)) > 1:
            auc = roc_auc_score(y_b, s_b)
            b_aucs.append(auc)
            b_strs.append(f"B{b}:{auc:.4f}")
        else:
            b_strs.append(f"B{b}:N/A (single class)")
            
    m_val = np.mean(b_aucs)
    sd_val = np.std(b_aucs, ddof=1)
    # 95% t-interval over valid blocks
    from scipy import stats
    df_t = len(b_aucs) - 1
    t_crit = stats.t.ppf(0.975, df_t)
    se = sd_val / np.sqrt(len(b_aucs))
    ci = [m_val - t_crit * se, m_val + t_crit * se]
    
    print(f"\n{name}:")
    print(f"  Per-block: {', '.join(b_strs)}")
    print(f"  Mean +/- SD: {m_val:.4f} +/- {sd_val:.4f}")
    print(f"  95% t-interval ({len(b_aucs)} blocks): [{ci[0]:.4f}, {ci[1]:.4f}]")

