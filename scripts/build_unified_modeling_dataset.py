import pathlib
"""
Assemble the Unified Modeling Dataset under Provisional Scope.
Strictly sourced from deterministic empirical data:
- HerpesDRG-DB (status == 'A', virus in [HSV1, HSV2])
- CMI 2025 Table S3 (Dähne et al. 2025)
Deduplicated on (Virus, Gene, Drug, Mutation).
Includes:
- Consensus phenotype & binary target (Resistant = 1, Sensitive = 0 at cited >= 3.0x cutoff)
- Flag for whether tuple is conflicted between sources
- Publication Year split: Pre-2016 (Train) vs Post-2016 (Prospective Test)
- 5-fold Position-Grouped cross-validation fold IDs
- Scope tier: Primary (UL23-ACV, UL30-ACV, UL30-FOS) vs Auxiliary / Exploratory
"""
import re
import json
import pandas as pd
import numpy as np
from sklearn.model_selection import GroupKFold

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
    if not isinstance(res_str, str): return 'Unknown'
    s = res_str.strip().lower()
    if 'resistant' in s and 'sensitive' not in s: return 'Resistant'
    if 'sensitive' in s and 'resistant' not in s: return 'Sensitive'
    if 'intermediate' in s: return 'Intermediate'
    if 'sensitive/resistant' in s: return 'Mixed/Borderline'
    return 'Unknown'

def normalize_herpesdrg_resistance(val):
    if pd.isna(val): return None
    s = str(val).strip()
    try:
        fc = float(s)
        return 'Resistant' if fc >= 3.0 else 'Sensitive'
    except ValueError:
        s_lower = s.lower()
        if 'resistant' in s_lower: return 'Resistant'
        if 'sensitive' in s_lower or 'susceptible' in s_lower or 'polymorphism' in s_lower: return 'Sensitive'
        return 'Unknown'

def extract_sources(row):
    sources = set()
    for col in ['pmids', 'doi', 'PMIDs', 'ref_doi', 'ref_link']:
        if col in row and pd.notna(row[col]):
            val = str(row[col]).strip()
            if val and val != 'nan':
                tokens = re.split(r'[\n\r\s,;]+', val)
                for t in tokens:
                    t = t.strip()
                    if t and t != 'nan':
                        sources.add(t)
    return sources

def main():
    hdrg_path = str(pathlib.Path(__file__).resolve().parent.parent / 'herpesdrg-db/herpesdrg-db.tsv')
    cmi_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/cmi_2025_all_sheets.csv')
    pub_years_path = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/publication_year_lookup.tsv')
    
    # 1. Load pub year lookup
    pub_df = pd.read_csv(pub_years_path, sep='\t')
    pub_map = dict(zip(pub_df['mutation_id'].astype(str), pub_df['paper_pub_year']))
    
    # 2. Load HerpesDRG
    hdrg = pd.read_csv(hdrg_path, sep='\t')
    hdrg = hdrg[(hdrg['status'] == 'A') & (hdrg['virus'].isin(['HSV1', 'HSV2']))].copy()
    hdrg['virus_norm'] = hdrg['virus'].map({'HSV1': 'HSV-1', 'HSV2': 'HSV-2'})
    
    hdrg_records = []
    drugs_list = ['Aciclovir', 'Foscarnet', 'Penciclovir', 'Brivudine', 'Cidofovir', 'Amenamevir', 'Pritelivir']
    for idx, r in hdrg.iterrows():
        v = r['virus_norm']
        g = r['gene']
        aa = str(r['aa_change']).strip()
        wt, pos, mut = parse_aa_change(aa)
        srcs = extract_sources(r)
        mut_id = str(r['mutation_id'])
        pyear = pub_map.get(mut_id, np.nan)
        
        for d in drugs_list:
            if pd.notna(r[d]):
                pheno = normalize_herpesdrg_resistance(r[d])
                hdrg_records.append({
                    'source_db': 'HerpesDRG',
                    'mutation_id': mut_id,
                    'virus': v,
                    'gene': g,
                    'drug': d,
                    'aa_change': aa,
                    'wt': wt,
                    'pos': pos,
                    'mut': mut,
                    'phenotype': pheno,
                    'raw_phenotype': r[d],
                    'pub_year': pyear,
                    'sources': srcs,
                    'strain_generation': r['strain_generation'],
                    'control_strain': r['control_strain']
                })
    hdrg_df = pd.DataFrame(hdrg_records)
    
    # 3. Load CMI 2025
    cmi = pd.read_csv(cmi_path)
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
        if g == 'Thymidine Kinase': gene_norm = 'UL23'
        elif g == 'DNA Polymerase': gene_norm = 'UL30'
        elif g in ['UL5', 'UL52']: gene_norm = g
        elif g == 'UL5/UL52': gene_norm = 'UL5/UL52'
        else: gene_norm = g
        
        drug = str(r['Drug']).strip()
        aa = str(r['AA change']).strip() if pd.notna(r['AA change']) else ''
        wt, pos_from_aa, mut = parse_aa_change(aa)
        res = normalize_cmi_resistance(r['Resistance'])
        pos_str = str(r['AA position']).strip()
        try: pos = int(pos_str)
        except: pos = pos_from_aa
        
        cmi_records.append({
            'source_db': 'CMI_2025',
            'mutation_id': f"CMI_{idx}",
            'virus': v,
            'gene': gene_norm,
            'drug': drug,
            'aa_change': aa,
            'wt': wt,
            'pos': pos,
            'mut': mut,
            'phenotype': res,
            'raw_phenotype': r['Resistance'],
            'pub_year': 2025, # CMI Review publication year
            'sources': extract_sources(r),
            'strain_generation': 'Review_Table',
            'control_strain': 'Reference'
        })
    cmi_df = pd.DataFrame(cmi_records)
    
    # 4. Detect conflicts on overlapping keys
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
            
    # Combine and deduplicate
    combined = pd.concat([hdrg_df, cmi_df], ignore_index=True)
    combined = combined[combined['pos'].notna()].copy()
    
    # Group by tuple_key to deduplicate
    dedup_rows = []
    for k, grp in combined.groupby('tuple_key'):
        v = grp['virus'].iloc[0]
        g = grp['gene'].iloc[0]
        d = grp['drug'].iloc[0]
        aa = grp['aa_change'].iloc[0]
        wt = grp['wt'].iloc[0]
        pos = int(grp['pos'].iloc[0])
        mut = grp['mut'].iloc[0]
        
        phenos = set(grp['phenotype'])
        is_conflicted = 1 if k in conflicted_tuples else 0
        
        # Consensus phenotype
        if is_conflicted:
            # Mark consensus as Ambiguous/Discordant
            consensus_pheno = 'Discordant'
            target_binary = np.nan
        elif 'Resistant' in phenos:
            consensus_pheno = 'Resistant'
            target_binary = 1
        elif 'Sensitive' in phenos:
            consensus_pheno = 'Sensitive'
            target_binary = 0
        else:
            consensus_pheno = 'Unknown'
            target_binary = np.nan
            
        # Earliest publication year among sources
        years = [float(y) for y in grp['pub_year'] if pd.notna(y)]
        min_year = min(years) if len(years) > 0 else np.nan
        
        # Temporal split (Pre-2016 vs Post-2016)
        if pd.notna(min_year):
            temporal_split = 'Pre-2016 (Train)' if min_year < 2016 else 'Post-2016 (Test)'
        else:
            temporal_split = 'Undetermined'
            
        all_srcs = set()
        for s_set in grp['sources']:
            all_srcs.update(s_set)
            
        # Primary scope eligibility
        is_primary = 1 if (g, d) in [('UL23', 'Aciclovir'), ('UL30', 'Aciclovir'), ('UL30', 'Foscarnet')] else 0
        
        dedup_rows.append({
            'tuple_key': k,
            'virus': v,
            'gene': g,
            'drug': d,
            'aa_change': aa,
            'wt': wt,
            'pos': pos,
            'mut': mut,
            'consensus_phenotype': consensus_pheno,
            'target_binary': target_binary,
            'is_conflicted': is_conflicted,
            'earliest_pub_year': min_year,
            'temporal_split': temporal_split,
            'is_primary_scope': is_primary,
            'num_sources': len(all_srcs),
            'sources': ';'.join(sorted(list(all_srcs)))
        })
        
    unified_df = pd.DataFrame(dedup_rows)
    
    # 5. Assign 5-fold Position-Grouped cross-validation fold IDs within each Gene-Drug pair
    unified_df['position_fold_id'] = -1
    gkf = GroupKFold(n_splits=5)
    
    for (g, d), grp_indices in unified_df.groupby(['gene', 'drug']).groups.items():
        sub_df = unified_df.loc[grp_indices]
        if sub_df['pos'].nunique() >= 5:
            for fold_i, (tr_idx, te_idx) in enumerate(gkf.split(sub_df, groups=sub_df['pos'])):
                actual_te_indices = grp_indices[te_idx]
                unified_df.loc[actual_te_indices, 'position_fold_id'] = fold_i
                
    out_tsv = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/unified_modeling_dataset.tsv')
    unified_df.to_csv(out_tsv, sep='\t', index=False)
    
    # Summary log
    log_lines = [
        "=== UNIFIED MODELING DATASET PROVENANCE & SUMMARY LOG ===",
        "Scope Status: Provisional Scope (Awaiting Human C1/C2 Tiering & PI Sign-Off)",
        f"Total Deduplicated Mutation-Drug Pairs: {len(unified_df)}",
        f"Primary Scope Rows (UL23-ACV, UL30-ACV, UL30-FOS): {len(unified_df[unified_df['is_primary_scope'] == 1])}",
        f"Discordant / Conflicted Tuples Flagged: {unified_df['is_conflicted'].sum()} ({unified_df['is_conflicted'].mean()*100:.2f}%)",
        "",
        "--- Primary Scope Summary ---"
    ]
    
    for pair in [('UL23', 'Aciclovir'), ('UL30', 'Aciclovir'), ('UL30', 'Foscarnet')]:
        p_sub = unified_df[(unified_df['gene'] == pair[0]) & (unified_df['drug'] == pair[1])]
        log_lines.append(f"Pair: {pair[0]}-{pair[1]}")
        log_lines.append(f"  Total Distinct Mutations: {len(p_sub)}")
        log_lines.append(f"  Distinct Positions: {p_sub['pos'].nunique()}")
        log_lines.append(f"  Resistant Mutations: {(p_sub['target_binary'] == 1).sum()}")
        log_lines.append(f"  Sensitive Mutations: {(p_sub['target_binary'] == 0).sum()}")
        log_lines.append(f"  Discordant (Excluded from clean train): {p_sub['is_conflicted'].sum()}")
        log_lines.append(f"  Pre-2016 Mutations: {(p_sub['temporal_split'] == 'Pre-2016 (Train)').sum()}")
        log_lines.append(f"  Post-2016 Mutations: {(p_sub['temporal_split'] == 'Post-2016 (Test)').sum()}")
        log_lines.append("")
        
    out_log = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/dataset_summary.log')
    with open(out_log, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))
        
    print("\n".join(log_lines))
    print(f"Saved dataset to: {out_tsv}")
    print(f"Saved log to: {out_log}")

if __name__ == '__main__':
    main()

