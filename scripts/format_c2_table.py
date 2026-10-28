import json
import pandas as pd

with open(r'd:\antiviral-resistance-hsv--1\logs\c2_split_evaluation_per_fold.json') as f:
    data = json.load(f)

rows = []
for pair, entries in data.items():
    for e in entries:
        rows.append({
            'Gene-Drug Pair': pair,
            'Regime': e['regime'],
            'Model': e['model'],
            'Mean Per-Fold AUROC': e['mean_per_fold_auroc'],
            'Std Per-Fold AUROC': e['std_per_fold_auroc'],
            'Pooled AUROC': e['pooled_auroc'],
            'CI_95_Low': e['ci_95'][0],
            'CI_95_High': e['ci_95'][1],
            'Per-Fold Folds': str(e['per_fold_values'])
        })

df = pd.DataFrame(rows)
out_tsv = r'd:\antiviral-resistance-hsv--1\logs\c2_full_models_regimes_table.tsv'
df.to_csv(out_tsv, sep='\t', index=False)
print("Saved full table to:", out_tsv)

# Print cleanly formatted markdown table for each pair
for pair in df['Gene-Drug Pair'].unique():
    print(f"\n### {pair} Benchmark Results (Mean Per-Fold +/- SD across 5 Folds)")
    sub = df[df['Gene-Drug Pair'] == pair]
    models = ['Logistic Regression', 'BLOSUM62 Radicalness', 'Position Prior', 'Physicochemical Random Forest']
    regimes = ['Regime 1: Random Split', 'Regime 2: Clade/Strain-Grouped Split', 'Regime 3: Position-Grouped Split', 'Regime 4: Both-Grouped Split']
    
    print(f"| Model | Random Split | Clade-Grouped | Position-Grouped | Both-Grouped |")
    print(f"| :--- | :--- | :--- | :--- | :--- |")
    for m in models:
        row_str = f"| **{m}** |"
        for r in regimes:
            match = sub[(sub['Model'] == m) & (sub['Regime'] == r)]
            if len(match) > 0:
                mean_v = match['Mean Per-Fold AUROC'].iloc[0]
                std_v = match['Std Per-Fold AUROC'].iloc[0]
                row_str += f" {mean_v:.4f} +/- {std_v:.4f} |"
            else:
                row_str += " N/A |"
        print(row_str)
