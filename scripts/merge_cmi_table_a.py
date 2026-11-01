import pathlib
"""
Merge CMI 2025 Table S3 with HerpesDRG-DB to compute unified Table A.
Detects overlaps, discordances/conflicts, and computes per-virus and combined counts.
"""
import re
import pandas as pd
import numpy as np

def parse_aa_change(aa_str):
    if not isinstance(aa_str, str):
        return None, None, None
    s = aa_str.strip()
    # Match standard single substitution like Y4D, D301N, A12P
    m = re.match(r'^([A-Z])(\d+)([A-Z*])$', s)
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    # Match stop codon like 8Stop
    m = re.match(r'^(\d+)(Stop|\*)$', s, re.IGNORECASE)
    if m:
        return 'X', int(m.group(1)), '*'
    # Match insertion / deletion / frameshift
    m = re.search(r'\d+', s)
    if m:
        return 'X', int(m.group()), 'X'
    return None, None, None

def normalize_cmi_resistance(res_str):
    if not isinstance(res_str, str):
        return 'Unknown'
    s = res_str.strip().lower()
    if 'resistant' in s and 'sensitive' not in s:
        return 'Resistant'
    if 'sensitive' in s and 'resistant' not in s:
        return 'Sensitive'
    if 'intermediate' in s:
        return 'Intermediate'
    if 'sensitive/resistant' in s:
        return 'Mixed/Borderline'
    return 'Unknown'

def normalize_herpesdrg_resistance(val):
    if pd.isna(val):
        return None
    s = str(val).strip()
    try:
        fc = float(s)
        return 'Resistant' if fc >= 3.0 else 'Sensitive'
    except ValueError:
        s_lower = s.lower()
        if 'resistant' in s_lower:
            return 'Resistant'
        if 'sensitive' in s_lower or 'susceptible' in s_lower or 'polymorphism' in s_lower:
            return 'Sensitive'
        return 'Unknown'

def main():
    hdrg_path = str(pathlib.Path(__file__).resolve().parent.parent / 'herpesdrg-db/herpesdrg-db.tsv')
    cmi_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/cmi_2025_all_sheets.csv')
    
    hdrg = pd.read_csv(hdrg_path, sep='\t')
    # Filter active HSV
    hdrg = hdrg[(hdrg['status'] == 'A') & (hdrg['virus'].isin(['HSV1', 'HSV2']))].copy()
    hdrg['virus_norm'] = hdrg['virus'].map({'HSV1': 'HSV-1', 'HSV2': 'HSV-2'})
    
    cmi = pd.read_csv(cmi_path)
    
    # Map CMI sheets to virus and gene
    sheet_virus_map = {
        'HSV-1 Pol': ('HSV-1', 'UL30'),
        'HSV-1 TK': ('HSV-1', 'UL23'),
        'HSV-2 Pol': ('HSV-2', 'UL30'),
        'HSV-2 TK': ('HSV-2', 'UL23'),
        'HSV-1 HPC': ('HSV-1', 'HPC'),
        'HSV-2 HPC': ('HSV-2', 'HPC')
    }
    
    cmi_records = []
    for idx, r in cmi.iterrows():
        sheet = r['Source_Sheet']
        v, g_default = sheet_virus_map.get(sheet, ('Unknown', 'Unknown'))
        g = r['Gene']
        if g == 'Thymidine Kinase':
            gene_norm = 'UL23'
        elif g == 'DNA Polymerase':
            gene_norm = 'UL30'
        elif g in ['UL5', 'UL52']:
            gene_norm = g
        elif g == 'UL5/UL52':
            gene_norm = 'UL5/UL52'
        else:
            gene_norm = g
            
        drug = str(r['Drug']).strip()
        aa_change = str(r['AA change']).strip() if pd.notna(r['AA change']) else ''
        res = normalize_cmi_resistance(r['Resistance'])
        pos_str = str(r['AA position']).strip()
        
        # parse position
        try:
            pos = int(pos_str)
        except:
            m = re.search(r'\d+', pos_str)
            pos = int(m.group()) if m else None
            
        cmi_records.append({
            'source': 'CMI_2025',
            'virus': v,
            'gene': gene_norm,
            'drug': drug,
            'aa_change': aa_change,
            'pos': pos,
            'phenotype': res,
            'raw_phenotype': r['Resistance'],
            'pmids': str(r['PMIDs']) if pd.notna(r['PMIDs']) else '',
            'culture_confirmed': r['resistance status confirmed in cell culture']
        })
    cmi_df = pd.DataFrame(cmi_records)
    
    # Process HerpesDRG into comparable long format
    hdrg_records = []
    drugs_list = ['Aciclovir', 'Foscarnet', 'Penciclovir', 'Brivudine', 'Cidofovir', 'Amenamevir', 'Pritelivir']
    for idx, r in hdrg.iterrows():
        v = r['virus_norm']
        g = r['gene']
        aa_change = str(r['aa_change']).strip()
        wt, pos, mut = parse_aa_change(aa_change)
        
        for d in drugs_list:
            if pd.notna(r[d]):
                pheno = normalize_herpesdrg_resistance(r[d])
                hdrg_records.append({
                    'source': 'HerpesDRG',
                    'virus': v,
                    'gene': g,
                    'drug': d,
                    'aa_change': aa_change,
                    'pos': pos,
                    'phenotype': pheno,
                    'raw_phenotype': r[d],
                    'pmids': str(r['ref_link']) if pd.notna(r['ref_link']) else '',
                    'doi': str(r['ref_doi']) if pd.notna(r['ref_doi']) else '',
                    'strain_gen': r['strain_generation']
                })
    hdrg_df = pd.DataFrame(hdrg_records)
    
    print(f"HerpesDRG active gene-drug entries: {len(hdrg_df)}")
    print(f"CMI 2025 active gene-drug entries: {len(cmi_df)}")
    
    # Find overlaps and conflicts
    # Key: (virus, gene, drug, aa_change)
    hdrg_df['key'] = hdrg_df['virus'] + '_' + hdrg_df['gene'] + '_' + hdrg_df['drug'] + '_' + hdrg_df['aa_change']
    cmi_df['key'] = cmi_df['virus'] + '_' + cmi_df['gene'] + '_' + cmi_df['drug'] + '_' + cmi_df['aa_change']
    
    common_keys = set(hdrg_df['key']).intersection(set(cmi_df['key']))
    print(f"Exact overlapping (Virus, Gene, Drug, Mutation) tuples: {len(common_keys)}")
    
    # Check concordance on overlapping keys
    conflicts = []
    concordant = 0
    hdrg_lookup = hdrg_df.groupby('key')['phenotype'].unique().to_dict()
    cmi_lookup = cmi_df.groupby('key')['phenotype'].unique().to_dict()
    
    for k in common_keys:
        h_p = hdrg_lookup[k]
        c_p = cmi_lookup[k]
        h_res = 'Resistant' in h_p
        c_res = 'Resistant' in c_p
        h_sen = 'Sensitive' in h_p
        c_sen = 'Sensitive' in c_p
        
        if (h_res and c_sen and not c_res) or (h_sen and c_res and not h_res):
            conflicts.append((k, h_p, c_p))
        else:
            concordant += 1
            
    print(f"Concordant tuples: {concordant}")
    print(f"Discordant / conflicting tuples: {len(conflicts)}")
    
    # Print examples of conflicts
    print("\n--- Sample Conflicts ---")
    for k, h_p, c_p in conflicts[:10]:
        print(f"  {k}: HerpesDRG={h_p} vs CMI={c_p}")
        
    # Build combined dataset
    combined = pd.concat([hdrg_df, cmi_df], ignore_index=True)
    
    # Now recount Table A across all Gene-Drug pairs:
    # 1. Combined (HSV-1 + HSV-2)
    # 2. Per-virus (HSV-1 vs HSV-2)
    pairs = [
        ('UL23', 'Aciclovir'),
        ('UL23', 'Foscarnet'),
        ('UL23', 'Penciclovir'),
        ('UL23', 'Brivudine'),
        ('UL30', 'Aciclovir'),
        ('UL30', 'Foscarnet'),
        ('UL30', 'Penciclovir'),
        ('UL30', 'Cidofovir'),
        ('UL5', 'Pritelivir'),
        ('UL5', 'Amenamevir'),
        ('UL52', 'Pritelivir'),
        ('UL52', 'Amenamevir'),
        ('UL29', 'Aciclovir')
    ]
    
    summary_rows = []
    for g, d in pairs:
        sub = combined[(combined['gene'] == g) & (combined['drug'] == d)].copy()
        sub = sub[sub['pos'].notna()]
        
        # Combined stats
        tot_rows = len(sub)
        hdrg_rows = len(sub[sub['source'] == 'HerpesDRG'])
        cmi_rows = len(sub[sub['source'] == 'CMI_2025'])
        
        distinct_pos = sub['pos'].nunique()
        res_pos = sub[sub['phenotype'] == 'Resistant']['pos'].nunique()
        sen_pos = sub[sub['phenotype'] == 'Sensitive']['pos'].nunique()
        
        # HSV-1 stats
        sub1 = sub[sub['virus'] == 'HSV-1']
        hsv1_pos = sub1['pos'].nunique()
        hsv1_res_pos = sub1[sub1['phenotype'] == 'Resistant']['pos'].nunique()
        
        # HSV-2 stats
        sub2 = sub[sub['virus'] == 'HSV-2']
        hsv2_pos = sub2['pos'].nunique()
        hsv2_res_pos = sub2[sub2['phenotype'] == 'Resistant']['pos'].nunique()
        
        # Audit Tier determination
        if res_pos >= 15:
            verdict = "PASS (>=15 res pos) - Primary"
        elif res_pos >= 5:
            verdict = "FAIL (<15 res pos) - Secondary / Pooled"
        else:
            verdict = "FAIL (<5 res pos) - Sparse / Relegated"
            
        summary_rows.append({
            'Gene': g,
            'Drug': d,
            'Total_Rows': tot_rows,
            'HerpesDRG_Rows': hdrg_rows,
            'CMI_Rows': cmi_rows,
            'Distinct_Positions': distinct_pos,
            'Resistant_Positions': res_pos,
            'Sensitive_Positions': sen_pos,
            'HSV1_Res_Positions': hsv1_res_pos,
            'HSV2_Res_Positions': hsv2_res_pos,
            'Verdict': verdict
        })
        
    sum_df = pd.DataFrame(summary_rows)
    print("\n=== RECOUNTED TABLE A (MERGED HERPESDRG + CMI 2025) ===")
    print(sum_df.to_string(index=False))
    
    # Save output
    sum_df.to_csv(str(pathlib.Path(__file__).resolve().parent.parent / 'logs/table_a_merged_recount.tsv'), sep='\t', index=False)
    
    # Write detailed conflicts log
    with open(str(pathlib.Path(__file__).resolve().parent.parent / 'logs/cmi_herpesdrg_conflicts.log'), 'w', encoding='utf-8') as f:
        f.write("=== CONFLICTS AND DISCORDANCE BETWEEN HERPESDRG-DB AND CMI 2025 TABLE S3 ===\n\n")
        f.write(f"Total overlapping tuples evaluated: {len(common_keys)}\n")
        f.write(f"Concordant tuples: {concordant} ({concordant/len(common_keys)*100:.2f}%)\n")
        f.write(f"Conflicting / discordant tuples: {len(conflicts)} ({len(conflicts)/len(common_keys)*100:.2f}%)\n\n")
        f.write("Detailed list of conflicts (HerpesDRG vs CMI 2025):\n")
        for k, h_p, c_p in conflicts:
            f.write(f"Tuple: {k}\n  HerpesDRG: {list(h_p)}\n  CMI 2025:  {list(c_p)}\n\n")

if __name__ == '__main__':
    main()

