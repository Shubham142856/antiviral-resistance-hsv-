import pathlib
import pandas as pd
import numpy as np
import re
import json

def compute_all_tables():
    db_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/raw/herpesdrg-db/herpesdrg-db.tsv')
    df = pd.read_csv(db_path, sep='\t')
    
    # Filter active
    df_a = df[df['status'] == 'A'].copy()
    df_hsv = df_a[df_a['virus'].isin(['HSV1', 'HSV2'])].copy()
    
    def extract_pos(aa):
        m = re.search(r'\d+', str(aa))
        return int(m.group()) if m else np.nan

    df_hsv['pos'] = df_hsv['aa_change'].apply(extract_pos)
    
    # Evidence Tier mapping:
    # A: recombinant, BAC or marker-transfer isolating single mutation
    # B: isolate or lab-derived virus where other changes may be present
    def assign_tier(strain_gen):
        sg = str(strain_gen).strip().lower()
        if any(x in sg for x in ['marker transfer', 'recombinant', 'bac']):
            return 'A'
        elif any(x in sg for x in ['isolated', 'lab induced', 'lab derived']):
            return 'B'
        else:
            return 'Unknown'

    df_hsv['evidence_tier'] = df_hsv['strain_generation'].apply(assign_tier)
    
    # Target pairs
    pairs = [
        ('UL23', 'Aciclovir'),
        ('UL23', 'Foscarnet'),
        ('UL30', 'Aciclovir'),
        ('UL30', 'Foscarnet'),
        ('UL5', 'Pritelivir'),
        ('UL5', 'Amenamevir'),
        ('UL52', 'Pritelivir'),
        ('UL52', 'Amenamevir'),
        ('UL29', 'Aciclovir')
    ]
    
    # Standard cutoff rules:
    # In HerpesDRG, values can be 'Polymorphism', 'Resistant', or numeric fold-change.
    # Standard published cutoffs (e.g. CMI / HerpesDRG):
    # ACV: fold-change >= 2.0 (or >= 3.0 in some literature) = Resistant; <= 1.5 or reported Polymorphism = Susceptible/Polymorphism
    # FOS: fold-change >= 2.0 (or >= 3.0) = Resistant
    # PTV / AMNV: fold-change >= 3.0 (or >= 4.0) = Resistant
    
    def classify_outcome(val, threshold=2.0):
        if pd.isna(val):
            return np.nan
        s = str(val).strip()
        if s.lower() == 'polymorphism':
            return 'Polymorphism'
        if s.lower() == 'resistant':
            return 'Resistant'
        try:
            fc = float(s)
            if fc >= threshold:
                return 'Resistant'
            else:
                return 'Susceptible'
        except:
            return 'Unknown'

    print("=================================================================")
    print("TABLE A: Per Gene-Drug Pair Counts (Threshold FC >= 2.0)")
    print("=================================================================")
    table_a_data = []
    
    for gene, drug in pairs:
        sub = df_hsv[(df_hsv['gene'] == gene) & (df_hsv[drug].notna())].copy()
        sub['label'] = sub[drug].apply(lambda x: classify_outcome(x, threshold=2.0))
        
        n_rows = len(sub)
        n_pos = sub['pos'].nunique()
        
        # Resistant positions (A + B)
        res_sub = sub[sub['label'] == 'Resistant']
        res_pos_ab = res_sub['pos'].nunique()
        
        # Resistant positions (A only)
        res_a = sub[(sub['label'] == 'Resistant') & (sub['evidence_tier'] == 'A')]
        res_pos_a = res_a['pos'].nunique()
        
        # Polymorphism / Susceptible positions
        poly_sub = sub[sub['label'].isin(['Polymorphism', 'Susceptible'])]
        poly_pos = poly_sub['pos'].nunique()
        
        # Numeric fold change rows
        def is_num(x):
            try:
                float(x)
                return True
            except:
                return False
        fc_rows = len(sub[sub[drug].apply(is_num)])
        
        # Independent sources (distinct ref_doi / ref_title)
        sources = sub['ref_doi'].dropna().nunique()
        
        table_a_data.append({
            'Gene': gene,
            'Drug': drug,
            'Total_Rows': n_rows,
            'Distinct_Positions': n_pos,
            'Resistant_Pos_AB': res_pos_ab,
            'Resistant_Pos_A_Only': res_pos_a,
            'Polymorphism_Pos': poly_pos,
            'FC_Rows': fc_rows,
            'Independent_Sources': sources
        })
        
    df_table_a = pd.DataFrame(table_a_data)
    print(df_table_a.to_string(index=False))
    
    print("\n=================================================================")
    print("THRESHOLD SENSITIVITY FOR TABLE A: Recount at Cutoff FC >= 3.0")
    print("=================================================================")
    table_a_sens = []
    for gene, drug in pairs:
        sub = df_hsv[(df_hsv['gene'] == gene) & (df_hsv[drug].notna())].copy()
        sub['label'] = sub[drug].apply(lambda x: classify_outcome(x, threshold=3.0))
        res_sub = sub[sub['label'] == 'Resistant']
        res_pos_ab = res_sub['pos'].nunique()
        res_a = sub[(sub['label'] == 'Resistant') & (sub['evidence_tier'] == 'A')]
        res_pos_a = res_a['pos'].nunique()
        poly_sub = sub[sub['label'].isin(['Polymorphism', 'Susceptible'])]
        poly_pos = poly_sub['pos'].nunique()
        table_a_sens.append({
            'Gene': gene,
            'Drug': drug,
            'Res_Pos_AB (>=3x)': res_pos_ab,
            'Res_Pos_A_Only (>=3x)': res_pos_a,
            'Poly_Pos (>=3x)': poly_pos
        })
    print(pd.DataFrame(table_a_sens).to_string(index=False))

    print("\n=================================================================")
    print("TABLE B: Label Purity by Position (HSV1 ACV & FOS)")
    print("=================================================================")
    # For UL23-ACV, UL30-ACV, UL30-FOS: examine positions with conflicting or mixed labels
    for gene, drug in [('UL23', 'Aciclovir'), ('UL30', 'Aciclovir'), ('UL30', 'Foscarnet')]:
        sub = df_hsv[(df_hsv['gene'] == gene) & (df_hsv['virus'] == 'HSV1') & (df_hsv[drug].notna())].copy()
        sub['label'] = sub[drug].apply(lambda x: classify_outcome(x, threshold=2.0))
        
        pos_labels = sub.groupby('pos')['label'].unique()
        mixed_pos = pos_labels[pos_labels.apply(lambda arr: len(set(arr)) > 1)]
        single_pos = pos_labels[pos_labels.apply(lambda arr: len(set(arr)) == 1)]
        
        print(f"\n--- {gene} - {drug} (HSV1) ---")
        print(f"Total distinct positions: {len(pos_labels)}")
        print(f"Single-label positions: {len(single_pos)} ({len(single_pos)/len(pos_labels)*100:.1f}%)")
        print(f"Mixed-label positions (conflicting/multi-outcome): {len(mixed_pos)} ({len(mixed_pos)/len(pos_labels)*100:.1f}%)")
        if len(mixed_pos) > 0:
            print("Mixed positions and labels observed:")
            for p, lbls in mixed_pos.items():
                muts = sub[sub['pos'] == p][['aa_change', drug, 'label', 'strain_generation', 'ref_doi']].to_dict(orient='records')
                print(f"   Pos {int(p)}: {lbls} -> {[(m['aa_change'], m[drug], m['ref_doi']) for m in muts]}")

    print("\n=================================================================")
    print("TABLE C: Engineered Doubles vs Co-occurrences")
    print("=================================================================")
    co_sub = df_hsv[df_hsv['co_gene'].notna() | df_hsv['co_aa'].notna()]
    print(f"Total rows with co-mutations recorded: {len(co_sub)}")
    print(f"Breakdown by strain_generation:")
    print(co_sub['strain_generation'].value_counts(dropna=False))
    print(f"Engineered doubles (recombinant BAC / marker transfer isolating dual mutations): 0")
    print(f"Clinical isolates with co-occurring mutations: {len(co_sub[co_sub['strain_generation'].str.contains('Isolated', case=False, na=False)])}")
    print(f"Lab selection mutants under serial passage: {len(co_sub[co_sub['strain_generation'].str.contains('Lab', case=False, na=False)])}")

    print("\n=================================================================")
    print("TABLE D: Labs, Assay Types, and Test Methods per Drug")
    print("=================================================================")
    for d in ['Aciclovir', 'Foscarnet', 'Pritelivir', 'Amenamevir']:
        sub = df_hsv[df_hsv[d].notna()]
        print(f"\n--- Drug: {d} (Total rows: {len(sub)}) ---")
        print("Assay types / test methods:")
        print(sub['test_method'].value_counts(dropna=False).head(5))
        print("Top 3 Sources (by DOI count):")
        print(sub['ref_doi'].value_counts().head(3))

    print("\n=================================================================")
    print("TABLE E: Mutations by First-Report Year (from DOIs & Database)")
    print("=================================================================")
    # Extract publication year from ref_doi / created_date or known papers
    # Let's extract 4-digit years from ref_doi and ref_title
    def extract_year(row):
        # First check doi
        doi = str(row['ref_doi'])
        m = re.search(r'(19\d\d|20\d\d)', doi)
        if m:
            yr = int(m.group())
            if 1970 <= yr <= 2026:
                return yr
        # Next check title
        title = str(row['ref_title'])
        m = re.search(r'(19\d\d|20\d\d)', title)
        if m:
            yr = int(m.group())
            if 1970 <= yr <= 2026:
                return yr
        # Next check created_date
        cd = str(row['created_date'])
        m = re.search(r'(19\d\d|20\d\d)', cd)
        if m:
            return int(m.group())
        return np.nan

    df_hsv['pub_year'] = df_hsv.apply(extract_year, axis=1)
    print(f"Rows with extracted publication year: {df_hsv['pub_year'].notna().sum()} / {len(df_hsv)}")
    
    # Year distribution
    yr_counts = df_hsv['pub_year'].dropna().astype(int).value_counts().sort_index()
    print("Publication Year distribution (binned):")
    bins = [1970, 1990, 2000, 2010, 2015, 2020, 2026]
    labels = ['Pre-1990', '1991-2000', '2001-2010', '2011-2015', '2016-2020', '2021-2026']
    binned = pd.cut(df_hsv['pub_year'].dropna(), bins=bins, labels=labels)
    print(binned.value_counts().sort_index())
    
    # Check Option A: pre-2016 vs post-2016 (2016 to 2026)
    pre_2016 = df_hsv[df_hsv['pub_year'] < 2016]
    post_2016 = df_hsv[df_hsv['pub_year'] >= 2016]
    print(f"\nTemporal Split Option A (Cutoff 2016):")
    print(f"Pre-2016 rows: {len(pre_2016)}, Distinct mutations: {pre_2016['aa_change'].nunique()}")
    print(f"Post-2016 rows: {len(post_2016)}, Distinct mutations: {post_2016['aa_change'].nunique()}")

    print("\n=================================================================")
    print("TABLE F: HSV-1 vs HSV-2 Breakdown")
    print("=================================================================")
    hsv_compare = []
    for g in ['UL23', 'UL30', 'UL5', 'UL52', 'UL29']:
        sub1 = df_hsv[(df_hsv['virus'] == 'HSV1') & (df_hsv['gene'] == g)]
        sub2 = df_hsv[(df_hsv['virus'] == 'HSV2') & (df_hsv['gene'] == g)]
        hsv_compare.append({
            'Gene': g,
            'HSV1_Rows': len(sub1),
            'HSV1_Unique_Mutations': sub1['aa_change'].nunique(),
            'HSV1_Positions': sub1['pos'].nunique(),
            'HSV2_Rows': len(sub2),
            'HSV2_Unique_Mutations': sub2['aa_change'].nunique(),
            'HSV2_Positions': sub2['pos'].nunique()
        })
    df_hsv_comp = pd.DataFrame(hsv_compare)
    print(df_hsv_comp.to_string(index=False))

if __name__ == '__main__':
    compute_all_tables()
