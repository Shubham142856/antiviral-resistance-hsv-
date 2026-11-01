import pathlib
import pandas as pd
import numpy as np
import hashlib
import re

def run_herpesdrg_audit():
    db_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/raw/herpesdrg-db/herpesdrg-db.tsv')
    
    with open(db_path, 'rb') as f:
        raw_bytes = f.read()
    raw_sha = hashlib.sha256(raw_bytes).hexdigest()
    lf_sha = hashlib.sha256(raw_bytes.replace(b'\r\n', b'\n')).hexdigest()
    
    print(f"=== HerpesDRG-DB File Hashes ===")
    print(f"Raw bytes: {len(raw_bytes)}, SHA256: {raw_sha}")
    print(f"LF normalized: {len(raw_bytes.replace(b'\r\n', b'\n'))}, SHA256: {lf_sha}")
    
    df = pd.read_csv(db_path, sep='\t')
    print(f"\nTotal rows in database: {len(df)}")
    print(f"Status value counts:\n{df['status'].value_counts(dropna=False)}")
    print(f"Virus value counts (all rows):\n{df['virus'].value_counts(dropna=False)}")
    
    df_a = df[df['status'] == 'A'].copy()
    print(f"\nActive rows (status == 'A'): {len(df_a)}")
    
    # Virus matching HSV1 and HSV2
    df_hsv = df_a[df_a['virus'].isin(['HSV1', 'HSV2'])].copy()
    print(f"Active HSV rows: {len(df_hsv)}")
    print(f"Active HSV1 rows: {len(df_hsv[df_hsv['virus'] == 'HSV1'])}")
    print(f"Active HSV2 rows: {len(df_hsv[df_hsv['virus'] == 'HSV2'])}")
    
    # Gene breakdown for active HSV
    print("\n=== Active HSV Rows by Virus and Gene ===")
    gene_summary = pd.crosstab(df_hsv['gene'], df_hsv['virus'], margins=True)
    print(gene_summary)
    
    # Drug columns
    drug_cols = ['Aciclovir', 'Foscarnet', 'Pritelivir', 'Amenamevir', 'Cidofovir', 
                 'Ganciclovir', 'Brincidofovir', 'Brivudine', 'Penciclovir']
    
    print("\n=== Active HSV Tested Rows per Drug ===")
    for d in drug_cols:
        tested = df_hsv[df_hsv[d].notna()]
        print(f"{d:15s}: {len(tested)} rows (HSV1: {len(tested[tested['virus']=='HSV1'])}, HSV2: {len(tested[tested['virus']=='HSV2'])})")
    
    # UL23 and Foscarnet investigation
    print("\n=== INVESTIGATION: UL23 and Foscarnet rows ===")
    ul23_fos = df_hsv[(df_hsv['gene'] == 'UL23') & (df_hsv['Foscarnet'].notna())]
    print(f"Total UL23 rows with Foscarnet value: {len(ul23_fos)}")
    for idx, r in ul23_fos.iterrows():
        print(f"ID: {r['mutation_id']}, Virus: {r['virus']}, AA: {r['aa_change']}, FOS_val: {r['Foscarnet']}, ACV_val: {r['Aciclovir']}, Strain: {r['strain_generation']}, Method: {r['test_method']}, Co-gene: {r['co_gene']}, Co-aa: {r['co_aa']}, DOI: {r['ref_doi']}, Title: {str(r['ref_title'])[:60]}..., Note: {r['note']}")

    # Table 5.1 replication
    print("\n=== Replication of Table 5.1 Baseline ===")
    def extract_pos(aa):
        m = re.search(r'\d+', str(aa))
        return int(m.group()) if m else np.nan

    df_hsv['pos'] = df_hsv['aa_change'].apply(extract_pos)
    
    target_pairs = [
        ('UL23', 'Aciclovir'),
        ('UL30', 'Aciclovir'),
        ('UL30', 'Foscarnet'),
        ('UL5', 'Pritelivir'),
        ('UL5', 'Amenamevir'),
        ('UL52', 'Pritelivir'),
        ('UL52', 'Amenamevir'),
        ('UL29', 'Aciclovir')
    ]
    
    for gene, drug in target_pairs:
        sub = df_hsv[(df_hsv['gene'] == gene) & (df_hsv[drug].notna())]
        n_rows = len(sub)
        distinct_pos = sub['pos'].nunique()
        
        def is_float(x):
            try:
                float(x)
                return True
            except:
                return False
                
        numeric_fc = sub[sub[drug].apply(is_float)]
        
        # Check values distribution
        vals = sub[drug].astype(str)
        # In HerpesDRG, what indicates resistance?
        # Many values are float fold-change, e.g. >= 2, >= 3, or strings
        print(f"\nPair ({gene}, {drug}): Rows={n_rows}, Distinct Pos={distinct_pos}, Numeric FC={len(numeric_fc)}")
        print(f"   Sample values / distribution:\n{vals.value_counts().head(8)}")

    # Let's inspect co-mutation rows
    print("\n=== Co-mutation rows in active HSV ===")
    co_rows = df_hsv[df_hsv['co_gene'].notna() | df_hsv['co_aa'].notna()]
    print(f"Total co-mutation rows: {len(co_rows)}")
    for idx, r in co_rows.iterrows():
        print(f"ID: {r['mutation_id']}, Virus: {r['virus']}, Gene: {r['gene']}, AA: {r['aa_change']}, Co-gene: {r['co_gene']}, Co-aa: {r['co_aa']}, Strain: {r['strain_generation']}, ACV: {r['Aciclovir']}, FOS: {r['Foscarnet']}, DOI: {r['ref_doi']}")

    # Let's inspect strain_generation values
    print("\n=== Strain generation / Evidence Tier breakdown ===")
    print(df_hsv['strain_generation'].value_counts(dropna=False))

    # Test method breakdown
    print("\n=== Test method / Assay breakdown ===")
    print(df_hsv['test_method'].value_counts(dropna=False))

    # Distinct mutations per gene (HSV1)
    print("\n=== Distinct unique mutations per gene (HSV1 active) ===")
    for g in df_hsv[df_hsv['virus']=='HSV1']['gene'].unique():
        subg = df_hsv[(df_hsv['virus']=='HSV1') & (df_hsv['gene']==g)]
        print(f"Gene {g}: {subg['aa_change'].nunique()} unique aa_changes, {len(subg)} rows")

if __name__ == '__main__':
    run_herpesdrg_audit()
