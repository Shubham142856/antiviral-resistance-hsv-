import pandas as pd
import numpy as np

def assign_contiguous_block(gene, pos):
    if gene == 'UL23':
        block = min(4, (pos - 1) // 75)
    else:
        block = min(4, (pos - 1) // 247)
    return int(block)

df = pd.read_csv(r'd:\antiviral-resistance-hsv--1\data\processed\primary_position_folds.tsv', sep='\t')
df['block'] = df.apply(lambda r: assign_contiguous_block(r['gene'], r['pos']), axis=1)

pairs = [
    ('UL23', 'Aciclovir', 'UL23-Aciclovir'),
    ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
    ('UL30', 'Foscarnet', 'UL30-Foscarnet')
]

print("=== CONTIGUOUS BLOCK COMPOSITION ACROSS TARGET PAIRS ===")
for g, d, name in pairs:
    sub = df[(df['gene'] == g) & (df['drug'] == d)]
    print(f"\nTarget Pair: {name} (Total rows: {len(sub)})")
    print(f"{'Block':<6} | {'Residue Range':<15} | {'Total Rows':<10} | {'Resistant':<10} | {'Sensitive':<10} | {'Res Ratio':<10}")
    print("-" * 75)
    for b in range(5):
        b_sub = sub[sub['block'] == b]
        n_res = (b_sub['target_binary'] == 1.0).sum()
        n_sen = (b_sub['target_binary'] == 0.0).sum()
        tot = len(b_sub)
        ratio = f"{n_res / tot:.1%}" if tot > 0 else "N/A"
        if g == 'UL23':
            r_str = f"{b*75 + 1} - {min(376, (b+1)*75)}"
        else:
            r_str = f"{b*247 + 1} - {min(1235, (b+1)*247)}"
        print(f"B{b:<5} | {r_str:<15} | {tot:<10} | {n_res:<10} | {n_sen:<10} | {ratio:<10}")

