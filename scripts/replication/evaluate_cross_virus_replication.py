"""
evaluate_cross_virus_replication.py
Evaluates rank-AUROC of BLOSUM62 radicalness across HCMV and VZV targets.
Anchored to pre-registration: docs/CROSS_VIRUS_REPLICATION_PREREGISTRATION.md
"""
import re
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

df = pd.read_csv('herpesdrg-db/herpesdrg-db.tsv', sep='\t')

# BLOSUM62 matrix
AA_LIST = list("ARNDCQEGHILKMFPSTWYV")
BLOSUM62 = {
    'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0},
    'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3},
    'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3},
    'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3},
    'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
    'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2},
    'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3},
    'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1},
    'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2},
    'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2},
    'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1},
    'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4}
}

def parse_aa(s):
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', str(s).strip())
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    return None

def is_res(val, threshold=2.0):
    s = str(val).strip()
    if s == 'Resistant': return 1
    if s in ['Polymorphism', 'Sensitive']: return 0
    try:
        return 1 if float(s) >= threshold else 0
    except:
        return None

def bootstrap_ci(y_true, scores, n_boot=2000, seed=42):
    rng = np.random.default_rng(seed)
    aucs = []
    n = len(y_true)
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        if len(np.unique(y_true[idx])) > 1:
            aucs.append(roc_auc_score(y_true[idx], scores[idx]))
    if len(aucs) == 0:
        return np.nan, np.nan
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))

pairs = [
    ('HCMV', 'UL97', 'Ganciclovir', 2.0, 'Kinase / Activator (Dispensable)'),
    ('HCMV', 'UL54', 'Ganciclovir', 2.0, 'DNA Polymerase (Essential)'),
    ('HCMV', 'UL54', 'Foscarnet', 2.0, 'DNA Polymerase (Essential)'),
    ('HCMV', 'UL54', 'Cidofovir', 2.0, 'DNA Polymerase (Essential)'),
    ('VZV', 'ORF36', 'Aciclovir', 3.0, 'Thymidine Kinase (Dispensable)'),
    ('VZV', 'ORF28', 'Aciclovir', 3.0, 'DNA Polymerase (Essential; Underpowered)')
]

results = {}

print("="*85)
print("CROSS-HERPESVIRUS REPLICATION OF RADICALNESS DISSOCIATION (HERPESDRG-DB)")
print("="*85)

for virus, gene, drug, thresh, role in pairs:
    sub = df[(df['virus'] == virus) & (df['gene'] == gene) & df[drug].notna()].copy()
    parsed = [parse_aa(x) for x in sub['aa_change']]
    sub['is_missense'] = [p is not None for p in parsed]
    sub['parsed'] = parsed
    miss = sub[sub['is_missense']].copy()
    
    labels = [is_res(x, thresh) for x in miss[drug]]
    valid_mask = [l is not None for l in labels]
    miss = miss[valid_mask].copy()
    miss['y'] = [l for l in labels if l is not None]
    
    # Deduplicate on unique (gene, aa_change)
    miss_dedup = miss.drop_duplicates(subset=['aa_change']).copy()
    
    # Compute BLOSUM62 radicalness (negative score: higher = more radical)
    rad_scores = []
    for idx_val, r in miss_dedup.iterrows():
        wt, pos, mut = r['parsed']
        # radicalness = - BLOSUM62(wt, mut)
        score = - BLOSUM62.get(wt, {}).get(mut, 0)
        rad_scores.append(score)
        
    miss_dedup['radicalness'] = rad_scores
    
    y = miss_dedup['y'].values
    scores = miss_dedup['radicalness'].values
    
    if len(np.unique(y)) > 1:
        auc = roc_auc_score(y, scores)
        ci_low, ci_high = bootstrap_ci(y, scores)
    else:
        auc, ci_low, ci_high = np.nan, np.nan, np.nan
        
    n_res = int((y == 1).sum())
    n_sens = int((y == 0).sum())
    
    results[f"{virus}_{gene}_{drug}"] = {
        'virus': virus,
        'gene': gene,
        'drug': drug,
        'role': role,
        'n_total': int(len(miss_dedup)),
        'n_res': n_res,
        'n_sens': n_sens,
        'blosum62_auroc': float(auc),
        'ci_95': [ci_low, ci_high]
    }
    
    print(f"\n{virus} {gene:<5} ({role}) - Drug: {drug}")
    print(f"  Sample: N={len(miss_dedup)} (Resistant={n_res}, Sensitive={n_sens})")
    print(f"  BLOSUM62 Radicalness AUROC: {auc:.4f} [95% CI: {ci_low:.4f}, {ci_high:.4f}]")

with open('logs/cross_virus_replication_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "="*85)
print("Results serialized to logs/cross_virus_replication_results.json")
print("="*85)
