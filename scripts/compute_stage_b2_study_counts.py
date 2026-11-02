"""
compute_stage_b2_study_counts.py
Stage B2: Print per-study counts using reference/PMID/DOI columns for each target
(N, resistant, sensitive per study) and the fraction of rows coming from single-class studies.
DO NOT DEFINE ANYTHING. STOP PROTOCOL ACTIVE.
"""
import re
import json
import numpy as np
import pandas as pd

# 1. Load primary benchmark cohort (929 rows)
primary_df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')

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

hsv1_herpes_dedup = hsv1_herpes.drop_duplicates(subset=['gene_std', 'pos', 'wt', 'mut']).copy()

# 3. Load CMI
cmi = pd.read_csv('data/processed/cmi_2025_all_sheets.csv')
cmi_hsv1 = cmi[cmi['Source_Sheet'].isin(['HSV-1 TK', 'HSV-1 Pol'])].copy()
cmi_hsv1['gene_std'] = cmi_hsv1['Gene'].astype(str).str.strip().replace({'Thymidine Kinase': 'UL23', 'DNA Polymerase': 'UL30'})
cmi_parsed = [parse_mut(x) for x in cmi_hsv1['AA change']]
cmi_hsv1['wt'] = [p[0] for p in cmi_parsed]
cmi_hsv1['pos'] = [p[1] for p in cmi_parsed]
cmi_hsv1['mut'] = [p[2] for p in cmi_parsed]
cmi_hsv1['drug_std'] = cmi_hsv1['Drug'].astype(str).str.strip()
cmi_hsv1_dedup = cmi_hsv1.drop_duplicates(subset=['gene_std', 'drug_std', 'pos', 'wt', 'mut']).copy()

# Merge
merged = pd.merge(
    primary_df,
    hsv1_herpes_dedup[['gene_std', 'pos', 'wt', 'mut', 'ref_doi', 'ref_title', 'ref_link']],
    left_on=['gene', 'pos', 'wt', 'mut'],
    right_on=['gene_std', 'pos', 'wt', 'mut'],
    how='left'
)

merged = pd.merge(
    merged,
    cmi_hsv1_dedup[['gene_std', 'drug_std', 'pos', 'wt', 'mut', 'PMIDs']],
    left_on=['gene', 'drug', 'pos', 'wt', 'mut'],
    right_on=['gene_std', 'drug_std', 'pos', 'wt', 'mut'],
    how='left'
)

# Unified study reference
def get_study_id(row):
    # Prefer DOI if available
    doi = str(row['ref_doi']).strip()
    if doi and doi not in ['nan', 'None', '']:
        return f"DOI:{doi}"
    pmid = str(row['PMIDs']).strip()
    if pmid and pmid not in ['nan', 'None', '']:
        return f"PMID:{pmid}"
    title = str(row['ref_title']).strip()
    if title and title not in ['nan', 'None', '']:
        return f"Title:{title[:40]}"
    return "Unassigned_or_Novel_Curation"

merged['study_id'] = [get_study_id(r) for _, r in merged.iterrows()]

targets = [
    ('UL30', 'Aciclovir', 'Primary Target'),
    ('UL30', 'Foscarnet', 'Secondary Target'),
    ('UL23', 'Aciclovir', 'Secondary Target')
]

results = {}

print("="*105)
print("STAGE B2: PER-STUDY COUNTS AND SINGLE-CLASS STUDY FRACTIONS")
print("="*105)

for g, d, desc in targets:
    t_name = f"{g}-{d}"
    sub = merged[(merged['gene'] == g) & (merged['drug'] == d)].copy()
    n_tot = len(sub)
    n_res = int(sub['target_binary'].sum())
    n_sens = n_tot - n_res
    
    # Group by study_id
    study_grp = sub.groupby('study_id')['target_binary'].agg(
        Total='count',
        Resistant='sum'
    ).reset_index()
    study_grp['Sensitive'] = study_grp['Total'] - study_grp['Resistant']
    study_grp['Pct_Resistant'] = study_grp['Resistant'] / study_grp['Total'] * 100
    
    # Classify study:
    # Single-class resistant: 100% resistant
    # Single-class sensitive: 100% sensitive (0% resistant)
    # Mixed: both classes present
    def classify_study(r):
        if r['Resistant'] == r['Total']:
            return 'All-Resistant'
        elif r['Sensitive'] == r['Total']:
            return 'All-Sensitive'
        else:
            return 'Mixed (Both Classes)'
            
    study_grp['Study_Type'] = [classify_study(r) for _, r in study_grp.iterrows()]
    study_grp = study_grp.sort_values(by='Total', ascending=False).reset_index(drop=True)
    
    n_studies = len(study_grp)
    n_single_class_studies = int((study_grp['Study_Type'] != 'Mixed (Both Classes)').sum())
    n_mixed_studies = n_studies - n_single_class_studies
    
    rows_single_class = int(study_grp[study_grp['Study_Type'] != 'Mixed (Both Classes)']['Total'].sum())
    frac_rows_single_class = rows_single_class / n_tot
    
    rows_all_res = int(study_grp[study_grp['Study_Type'] == 'All-Resistant']['Total'].sum())
    rows_all_sens = int(study_grp[study_grp['Study_Type'] == 'All-Sensitive']['Total'].sum())
    rows_mixed = int(study_grp[study_grp['Study_Type'] == 'Mixed (Both Classes)']['Total'].sum())
    
    print(f"\n################################################################################")
    print(f"TARGET: {t_name} ({desc}) | N = {n_tot} (Res={n_res}, Sens={n_sens})")
    print(f"Total Studies: {n_studies} | Single-Class Studies: {n_single_class_studies} ({n_single_class_studies/n_studies*100:.1f}%) | Mixed Studies: {n_mixed_studies}")
    print(f"Fraction of rows from Single-Class Studies: {rows_single_class}/{n_tot} = {frac_rows_single_class*100:.2f}%")
    print(f"  - Rows from All-Sensitive Studies: {rows_all_sens} ({rows_all_sens/n_tot*100:.1f}%)")
    print(f"  - Rows from All-Resistant Studies: {rows_all_res} ({rows_all_res/n_tot*100:.1f}%)")
    print(f"  - Rows from Mixed (Both-Class) Studies: {rows_mixed} ({rows_mixed/n_tot*100:.1f}%)")
    print("-" * 105)
    print(f"{'Study Identifier / DOI / PMID':<55} | {'Total':<6} | {'Res':<5} | {'Sens':<5} | {'Study Type'}")
    print("-" * 105)
    
    studies_list = []
    for _, r in study_grp.iterrows():
        print(f"{str(r['study_id'])[:55]:<55} | {r['Total']:<6} | {r['Resistant']:<5} | {r['Sensitive']:<5} | {r['Study_Type']}")
        studies_list.append({
            'study_id': str(r['study_id']),
            'total': int(r['Total']),
            'resistant': int(r['Resistant']),
            'sensitive': int(r['Sensitive']),
            'pct_resistant': float(r['Pct_Resistant']),
            'study_type': r['Study_Type']
        })
        
    results[t_name] = {
        'total_rows': n_tot,
        'resistant_rows': n_res,
        'sensitive_rows': n_sens,
        'total_studies': n_studies,
        'single_class_studies_count': n_single_class_studies,
        'mixed_studies_count': n_mixed_studies,
        'fraction_rows_from_single_class_studies': frac_rows_single_class,
        'rows_from_all_sensitive_studies': rows_all_sens,
        'rows_from_all_resistant_studies': rows_all_res,
        'rows_from_mixed_studies': rows_mixed,
        'studies': studies_list
    }

with open('logs/stage_b2_study_counts.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "="*105)
print("Saved Stage B2 breakdown to logs/stage_b2_study_counts.json")
print("STAGE B2 COMPLETE. STOPPING AS DIRECTED.")
print("="*105)
