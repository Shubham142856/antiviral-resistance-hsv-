import pathlib
"""
Deduplicated Table A with Conflict Handling (With vs Without Conflicted Tuples).
Deduplicates on (Virus, Gene, Drug, Mutation).
Reports distinct mutations, distinct positions, and distinct source publications.
Compares counts including vs strictly excluding the 65 discordant tuples.
"""
import re
import pandas as pd
import numpy as np

def parse_aa_change(aa_str):
    if not isinstance(aa_str, str):
        return None, None, None
    s = aa_str.strip()
    m = re.match(r'^([A-Z])(\d+)([A-Z*])$', s)
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    m = re.match(r'^(\d+)(Stop|\*)$', s, re.IGNORECASE)
    if m:
        return 'X', int(m.group(1)), '*'
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

def extract_sources(row):
    sources = set()
    for col in ['pmids', 'doi', 'PMIDs', 'ref_doi', 'ref_link']:
        if col in row and pd.notna(row[col]):
            val = str(row[col]).strip()
            if val and val != 'nan':
                # split by newline, space, or comma
                tokens = re.split(r'[\n\r\s,;]+', val)
                for t in tokens:
                    t = t.strip()
                    if t and t != 'nan':
                        sources.add(t)
    return sources

def main():
    hdrg_path = str(pathlib.Path(__file__).resolve().parent.parent / 'herpesdrg-db/herpesdrg-db.tsv')
    cmi_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/cmi_2025_all_sheets.csv')
    
    hdrg = pd.read_csv(hdrg_path, sep='\t')
    hdrg = hdrg[(hdrg['status'] == 'A') & (hdrg['virus'].isin(['HSV1', 'HSV2']))].copy()
    hdrg['virus_norm'] = hdrg['virus'].map({'HSV1': 'HSV-1', 'HSV2': 'HSV-2'})
    
    cmi = pd.read_csv(cmi_path)
    sheet_virus_map = {
        'HSV-1 Pol': ('HSV-1', 'UL30'),
        'HSV-1 TK': ('HSV-1', 'UL23'),
        'HSV-2 Pol': ('HSV-2', 'UL30'),
        'HSV-2 TK': ('HSV-2', 'UL23'),
        'HSV-1 HPC': ('HSV-1', 'HPC'),
        'HSV-2 HPC': ('HSV-2', 'HPC')
    }
    
    # 1. Standardize CMI records
    cmi_rows = []
    for idx, r in cmi.iterrows():
        sheet = r['Source_Sheet']
        v, g_default = sheet_virus_map.get(sheet, ('Unknown', 'Unknown'))
        g = r['Gene']
        if g == 'Thymidine Kinase': gene_norm = 'UL23'
        elif g == 'DNA Polymerase': gene_norm = 'UL30'
        elif g in ['UL5', 'UL52']: gene_norm = g
        elif g == 'UL5/UL52': gene_norm = 'UL5/UL52'
        else: gene_norm = g
        
        drug = str(r['Drug']).strip()
        aa_change = str(r['AA change']).strip() if pd.notna(r['AA change']) else ''
        res = normalize_cmi_resistance(r['Resistance'])
        pos_str = str(r['AA position']).strip()
        try: pos = int(pos_str)
        except:
            m = re.search(r'\d+', pos_str)
            pos = int(m.group()) if m else None
            
        cmi_rows.append({
            'source_db': 'CMI_2025',
            'virus': v,
            'gene': gene_norm,
            'drug': drug,
            'aa_change': aa_change,
            'pos': pos,
            'phenotype': res,
            'raw_phenotype': r['Resistance'],
            'sources': extract_sources(r)
        })
    cmi_df = pd.DataFrame(cmi_rows)
    
    # 2. Standardize HerpesDRG records
    hdrg_rows = []
    drugs_list = ['Aciclovir', 'Foscarnet', 'Penciclovir', 'Brivudine', 'Cidofovir', 'Amenamevir', 'Pritelivir']
    for idx, r in hdrg.iterrows():
        v = r['virus_norm']
        g = r['gene']
        aa_change = str(r['aa_change']).strip()
        wt, pos, mut = parse_aa_change(aa_change)
        srcs = extract_sources(r)
        
        for d in drugs_list:
            if pd.notna(r[d]):
                pheno = normalize_herpesdrg_resistance(r[d])
                hdrg_rows.append({
                    'source_db': 'HerpesDRG',
                    'virus': v,
                    'gene': g,
                    'drug': d,
                    'aa_change': aa_change,
                    'pos': pos,
                    'phenotype': pheno,
                    'raw_phenotype': r[d],
                    'sources': srcs
                })
    hdrg_df = pd.DataFrame(hdrg_rows)
    
    # 3. Identify overlapping tuples and detect conflicts
    hdrg_df['tuple_key'] = hdrg_df['virus'] + '|' + hdrg_df['gene'] + '|' + hdrg_df['drug'] + '|' + hdrg_df['aa_change']
    cmi_df['tuple_key'] = cmi_df['virus'] + '|' + cmi_df['gene'] + '|' + cmi_df['drug'] + '|' + cmi_df['aa_change']
    
    common_tuples = set(hdrg_df['tuple_key']).intersection(set(cmi_df['tuple_key']))
    
    conflicted_tuples = set()
    for k in common_tuples:
        h_phenos = set(hdrg_df[hdrg_df['tuple_key'] == k]['phenotype'])
        c_phenos = set(cmi_df[cmi_df['tuple_key'] == k]['phenotype'])
        h_res = 'Resistant' in h_phenos
        c_res = 'Resistant' in c_phenos
        h_sen = 'Sensitive' in h_phenos
        c_sen = 'Sensitive' in c_phenos
        if (h_res and c_sen and not c_res) or (h_sen and c_res and not h_res):
            conflicted_tuples.add(k)
            
    print(f"Total overlapping tuples: {len(common_tuples)}")
    print(f"Total conflicted / discordant tuples: {len(conflicted_tuples)}")
    
    # 4. Deduplicate into unique mutation tuples
    # For each tuple, combine sources; assign consensus phenotype
    all_combined = pd.concat([hdrg_df, cmi_df], ignore_index=True)
    all_combined = all_combined[all_combined['pos'].notna()].copy()
    
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
    
    results = []
    for g, d in pairs:
        sub = all_combined[(all_combined['gene'] == g) & (all_combined['drug'] == d)].copy()
        
        # Breakdown: All vs Clean (excluding conflicted)
        for mode in ['All_Inclusive', 'Strict_Excluding_Conflicts']:
            if mode == 'Strict_Excluding_Conflicts':
                sub_mode = sub[~sub['tuple_key'].isin(conflicted_tuples)].copy()
            else:
                sub_mode = sub.copy()
                
            # Deduplicate by tuple_key
            # Group by tuple_key: aggregate phenotypes and sources
            def agg_tuple(group):
                phenos = set(group['phenotype'])
                res = 'Resistant' if 'Resistant' in phenos else ('Sensitive' if 'Sensitive' in phenos else 'Unknown')
                all_srcs = set()
                for s_set in group['sources']:
                    all_srcs.update(s_set)
                return pd.Series({
                    'virus': group['virus'].iloc[0],
                    'pos': group['pos'].iloc[0],
                    'consensus_phenotype': res,
                    'num_sources': len(all_srcs),
                    'sources_list': ';'.join(sorted(list(all_srcs)))
                })
                
            dedup = sub_mode.groupby('tuple_key').apply(agg_tuple).reset_index()
            
            # Counts
            total_distinct_muts = len(dedup)
            total_distinct_positions = dedup['pos'].nunique()
            
            res_muts = len(dedup[dedup['consensus_phenotype'] == 'Resistant'])
            res_pos = dedup[dedup['consensus_phenotype'] == 'Resistant']['pos'].nunique()
            
            sen_muts = len(dedup[dedup['consensus_phenotype'] == 'Sensitive'])
            sen_pos = dedup[dedup['consensus_phenotype'] == 'Sensitive']['pos'].nunique()
            
            # Distinct sources
            all_pair_sources = set()
            for s_str in dedup['sources_list']:
                if s_str:
                    all_pair_sources.update(s_str.split(';'))
            distinct_sources_count = len(all_pair_sources)
            
            # Per-virus resistant positions
            dedup_hsv1 = dedup[dedup['virus'] == 'HSV-1']
            hsv1_res_pos = dedup_hsv1[dedup_hsv1['consensus_phenotype'] == 'Resistant']['pos'].nunique()
            hsv1_res_muts = len(dedup_hsv1[dedup_hsv1['consensus_phenotype'] == 'Resistant'])
            
            dedup_hsv2 = dedup[dedup['virus'] == 'HSV-2']
            hsv2_res_pos = dedup_hsv2[dedup_hsv2['consensus_phenotype'] == 'Resistant']['pos'].nunique()
            hsv2_res_muts = len(dedup_hsv2[dedup_hsv2['consensus_phenotype'] == 'Resistant'])
            
            if res_pos >= 15:
                tier = "PASS (>=15 res pos) - Primary"
            elif res_pos >= 5:
                tier = "FAIL (<15 res pos) - Secondary / Pooled"
            else:
                tier = "FAIL (<5 res pos) - Relegated / Underpowered"
                
            results.append({
                'Gene': g,
                'Drug': d,
                'Evaluation_Mode': mode,
                'Distinct_Mutations': total_distinct_muts,
                'Distinct_Positions': total_distinct_positions,
                'Resistant_Mutations': res_muts,
                'Resistant_Positions': res_pos,
                'Sensitive_Positions': sen_pos,
                'HSV1_Res_Positions': hsv1_res_pos,
                'HSV2_Res_Positions': hsv2_res_pos,
                'Distinct_Literature_Sources': distinct_sources_count,
                'Audit_Tier': tier
            })
            
    res_df = pd.DataFrame(results)
    print("\n=== DEDUPLICATED TABLE A: WITH VS WITHOUT CONFLICTED TUPLES ===")
    print(res_df.to_string(index=False))
    
    out_tsv = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/table_a_deduplicated.tsv')
    res_df.to_csv(out_tsv, sep='\t', index=False)
    print(f"\nSaved deduplicated Table A to: {out_tsv}")

if __name__ == '__main__':
    main()

