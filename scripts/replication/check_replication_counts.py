"""
check_replication_counts.py
Parses missense mutations and checks class counts for HCMV and VZV in HerpesDRG-DB.
"""
import re
import pandas as pd

df = pd.read_csv('herpesdrg-db/herpesdrg-db.tsv', sep='\t')

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

pairs = [
    ('HCMV', 'UL97', 'Ganciclovir', 2.0),
    ('HCMV', 'UL54', 'Ganciclovir', 2.0),
    ('HCMV', 'UL54', 'Foscarnet', 2.0),
    ('HCMV', 'UL54', 'Cidofovir', 2.0),
    ('VZV', 'ORF36', 'Aciclovir', 3.0),
    ('VZV', 'ORF28', 'Aciclovir', 3.0),
    ('VZV', 'ORF28', 'Brivudine', 3.0)
]

print("="*75)
print("HERPESDRG-DB REPLICATION COHORT COUNTS")
print("="*75)

for virus, gene, drug, thresh in pairs:
    sub = df[(df['virus'] == virus) & (df['gene'] == gene) & df[drug].notna()].copy()
    parsed = [parse_aa(x) for x in sub['aa_change']]
    sub['is_missense'] = [p is not None for p in parsed]
    miss = sub[sub['is_missense']].copy()
    
    labels = [is_res(x, thresh) for x in miss[drug]]
    valid_mask = [l is not None for l in labels]
    miss = miss[valid_mask].copy()
    miss['y'] = [l for l in labels if l is not None]
    
    # Deduplicate on unique (gene, pos, mut)
    miss_dedup = miss.drop_duplicates(subset=['aa_change'])
    
    n_res = int((miss_dedup['y'] == 1).sum())
    n_sens = int((miss_dedup['y'] == 0).sum())
    n_tot = len(miss_dedup)
    
    print(f"{virus} {gene:<5} - {drug:<12} (Cutoff >={thresh}x): Total={n_tot:<3} | Res={n_res:<3} | Sens={n_sens:<3} | Evaluable: {n_res >= 10 and n_sens >= 10}")
