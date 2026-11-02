"""
inspect_all_evidence_columns.py
STAGE B: Print unique values and counts of every candidate evidence column
in HerpesDRG-DB and CMI 2025 (all HSV-1 sheets) per target and by resistant vs sensitive.
DO NOT DEFINE ANY TIER MAPPING.
"""
import re
import json
import numpy as np
import pandas as pd

# 1. Load primary benchmark cohort (929 rows)
primary_df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
print(f"Loaded primary benchmark cohort: N={len(primary_df)} rows")

# 2. Load HerpesDRG
herpes = pd.read_csv('herpesdrg-db/herpesdrg-db.tsv', sep='\t')
hsv1_herpes = herpes[herpes['virus'] == 'HSV1'].copy()

def parse_mut(s):
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', str(s).strip())
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    return None, None, None

parsed = [parse_mut(x) for x in hsv1_herpes['aa_change']]
hsv1_herpes['wt'] = [p[0] for p in parsed]
hsv1_herpes['pos'] = [p[1] for p in parsed]
hsv1_herpes['mut'] = [p[2] for p in parsed]
hsv1_herpes['gene_std'] = hsv1_herpes['gene'].astype(str).str.strip().replace({'TK': 'UL23', 'POL': 'UL30', 'DPOL': 'UL30'})

# Strip whitespace from string columns in herpes
for col in ['strain_generation', 'test_method', 'control_strain']:
    hsv1_herpes[col] = hsv1_herpes[col].astype(str).str.strip()
    hsv1_herpes.loc[hsv1_herpes[col].isin(['nan', 'None', '']), col] = np.nan

hdrg_cols = ['strain_generation', 'test_method', 'control_strain', 'Aciclovir', 'Foscarnet', 'ref_doi']
hsv1_herpes_dedup = hsv1_herpes.drop_duplicates(subset=['gene_std', 'pos', 'wt', 'mut']).copy()

# 3. Load CMI 2025 All Sheets (HSV-1 TK and HSV-1 Pol)
cmi = pd.read_csv('data/processed/cmi_2025_all_sheets.csv')
cmi_hsv1 = cmi[cmi['Source_Sheet'].isin(['HSV-1 TK', 'HSV-1 Pol'])].copy()
cmi_hsv1['gene_std'] = cmi_hsv1['Gene'].astype(str).str.strip().replace({
    'Thymidine Kinase': 'UL23',
    'DNA Polymerase': 'UL30'
})
cmi_parsed = [parse_mut(x) for x in cmi_hsv1['AA change']]
cmi_hsv1['wt'] = [p[0] for p in cmi_parsed]
cmi_hsv1['pos'] = [p[1] for p in cmi_parsed]
cmi_hsv1['mut'] = [p[2] for p in cmi_parsed]
cmi_hsv1['drug_std'] = cmi_hsv1['Drug'].astype(str).str.strip()

for col in ['resistance status confirmed in cell culture', 'clinical susceptibility information available', 'All genes sequenced that may confer resistance?', 'Resistance']:
    cmi_hsv1[col] = cmi_hsv1[col].astype(str).str.strip()
    cmi_hsv1.loc[cmi_hsv1[col].isin(['nan', 'None', '']), col] = np.nan

cmi_cols = [
    'resistance status confirmed in cell culture',
    'clinical susceptibility information available',
    'All genes sequenced that may confer resistance?',
    'Resistance'
]
cmi_hsv1_dedup = cmi_hsv1.drop_duplicates(subset=['gene_std', 'drug_std', 'pos', 'wt', 'mut']).copy()

# Merge onto primary_df
merged = pd.merge(
    primary_df,
    hsv1_herpes_dedup[['gene_std', 'pos', 'wt', 'mut'] + hdrg_cols],
    left_on=['gene', 'pos', 'wt', 'mut'],
    right_on=['gene_std', 'pos', 'wt', 'mut'],
    how='left'
)

merged = pd.merge(
    merged,
    cmi_hsv1_dedup[['gene_std', 'drug_std', 'pos', 'wt', 'mut'] + cmi_cols],
    left_on=['gene', 'drug', 'pos', 'wt', 'mut'],
    right_on=['gene_std', 'drug_std', 'pos', 'wt', 'mut'],
    how='left'
)

# Helper to check if drug value in HerpesDRG is numeric
def get_hdrg_drug_type(row):
    val = row['Aciclovir'] if row['drug'] == 'Aciclovir' else row['Foscarnet']
    if pd.isna(val) or str(val).strip() in ['nan', 'None', '']:
        return '<Missing in HerpesDRG>'
    s = str(val).replace('>', '').replace('<', '').strip()
    try:
        float(s)
        return 'Numeric Fold-Change (e.g. 2.5, 5.0)'
    except:
        return f"Categorical String: '{val.strip()}'"

merged['herpesdrg_phenotype_format'] = [get_hdrg_drug_type(r) for _, r in merged.iterrows()]

candidate_cols = [
    ('strain_generation', 'HerpesDRG: Strain Generation / Origin'),
    ('test_method', 'HerpesDRG: Test Method / Phenotypic Assay'),
    ('control_strain', 'HerpesDRG: Reference Control Strain'),
    ('herpesdrg_phenotype_format', 'HerpesDRG: Drug Column Format (Numeric Fold-Change vs Categorical String)'),
    ('resistance status confirmed in cell culture', 'CMI 2025: Resistance Confirmed in Cell Culture'),
    ('clinical susceptibility information available', 'CMI 2025: Clinical Susceptibility Info Available'),
    ('All genes sequenced that may confer resistance?', 'CMI 2025: All Candidate Resistance Genes Sequenced'),
    ('Resistance', 'CMI 2025: Resistance Phenotype Label')
]

targets = [
    ('UL30', 'Aciclovir', 'Primary Target'),
    ('UL30', 'Foscarnet', 'Secondary Target'),
    ('UL23', 'Aciclovir', 'Secondary Target')
]

print("\n" + "="*115)
print("STAGE B: RAW CANDIDATE EVIDENCE COLUMNS (VALUES & COUNTS BY TARGET AND RESISTANCE LABEL)")
print("RULE ENFORCED: NO TIER MAPPING HAS BEEN APPLIED OR DEFINED. STOP PROTOCOL ACTIVE.")
print("="*115)

all_stats = {}

for g, d, desc in targets:
    t_name = f"{g}-{d}"
    sub = merged[(merged['gene'] == g) & (merged['drug'] == d)].copy()
    n_tot = len(sub)
    n_res = int(sub['target_binary'].sum())
    n_sens = n_tot - n_res
    print(f"\n" + "#"*115)
    print(f"TARGET: {t_name} ({desc}) | Total Records = {n_tot} (Resistant = {n_res}, Sensitive = {n_sens})")
    print("#"*115)
    
    t_data = {'total': n_tot, 'resistant': n_res, 'sensitive': n_sens, 'columns': {}}
    
    for col_name, col_desc in candidate_cols:
        print(f"\n--- [Column: {col_name}] ({col_desc}) ---")
        sub_col = sub[col_name].fillna('<MISSING / NOT PRESENT>')
        
        ct = pd.crosstab(sub_col, sub['target_binary'])
        if 0 not in ct.columns: ct[0] = 0
        if 1 not in ct.columns: ct[1] = 0
        ct = ct.rename(columns={0: 'Sensitive', 1: 'Resistant'})
        ct['Total'] = ct['Sensitive'] + ct['Resistant']
        ct = ct.sort_values(by='Total', ascending=False)
        
        print(f"{'Value':<48} | {'Total':<8} | {'Resistant':<10} | {'Sensitive':<10} | {'% Res':<8}")
        print("-" * 92)
        c_dict = {}
        for val, row in ct.iterrows():
            pct_res = (row['Resistant'] / row['Total'] * 100) if row['Total'] > 0 else 0
            print(f"{str(val):<48} | {row['Total']:<8} | {row['Resistant']:<10} | {row['Sensitive']:<10} | {pct_res:5.1f}%")
            c_dict[str(val)] = {
                'total': int(row['Total']),
                'resistant': int(row['Resistant']),
                'sensitive': int(row['Sensitive']),
                'pct_resistant': float(pct_res)
            }
        t_data['columns'][col_name] = c_dict
        
    all_stats[t_name] = t_data

with open('logs/stage_b_candidate_evidence_columns.json', 'w') as f:
    json.dump(all_stats, f, indent=2)

print("\n" + "="*115)
print("Saved complete JSON breakdown to logs/stage_b_candidate_evidence_columns.json")
print("STAGE B COMPLETE. STOPPING HERE FOR USER TIER MAPPING DEFINITION.")
print("="*115)
