import pathlib
"""
Build Reproducible Folds, Multi-Mutation Exclusions, and Dataset Enhancements:
1. Enhances data/processed/unified_modeling_dataset.tsv:
   - Adds 'is_single_site' (1 or 0)
   - Adds 'mutation_category' ('single_substitution', 'nonsense_stop', 'frameshift', 'indel', 'co_occurrence_combo', 'other')
2. Assembles data/processed/co_occurrence_entries.tsv:
   - Co-mutations from HerpesDRG (co_gene, co_aa)
   - Combination entries from CMI Table S3 ('/', '+', ';', multi-gene)
3. Exports data/processed/primary_position_folds.tsv:
   - Primary scope (HSV-1 UL23-ACV, UL30-ACV, UL30-FOS)
   - Single-site substitutions only (is_single_site == 1)
   - Deterministic 5-fold Position-Grouped cross-validation (GroupKFold, seed=42)
   - Hard assertion: zero position leakage across train/test splits in every fold
4. Computes SHA-256 hash of primary_position_folds.tsv and records in logs/fold_assignments_hash.log
"""
import re
import hashlib
import pandas as pd
import numpy as np
from sklearn.model_selection import GroupKFold

CANONICAL_AAS = set('ACDEFGHIKLMNPQRSTVWY')

def categorize_mutation(aa_str):
    s = str(aa_str).strip()
    # Check combo / co-occurrence
    if '/' in s or '+' in s or ';' in s:
        return 'co_occurrence_combo', 0
    # Check pure single substitution
    m = re.match(r'^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$', s)
    if m:
        wt, pos, mut = m.groups()
        if wt != mut:
            return 'single_substitution', 1
        else:
            return 'synonymous', 0
    # Check stop / nonsense
    if '*' in s or 'stop' in s.lower():
        return 'nonsense_stop', 0
    # Check frameshift
    if 'fs' in s.lower() or 'frameshift' in s.lower():
        return 'frameshift', 0
    # Check indel
    if 'del' in s.lower() or 'ins' in s.lower():
        return 'indel', 0
    return 'other', 0

def main():
    ds_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/unified_modeling_dataset.tsv')
    df = pd.read_csv(ds_path, sep='\t')
    
    # 1. Add categories and single site flag
    cats = []
    singles = []
    for idx, r in df.iterrows():
        cat, is_s = categorize_mutation(r['aa_change'])
        cats.append(cat)
        singles.append(is_s)
        
    df['mutation_category'] = cats
    df['is_single_site'] = singles
    
    # Overwrite unified modeling dataset
    df.to_csv(ds_path, sep='\t', index=False)
    print(f"Updated unified dataset at {ds_path} with is_single_site and mutation_category.")
    print("Category counts:")
    print(df['mutation_category'].value_counts())
    
    # 2. Extract co-occurrence entries from HerpesDRG and CMI
    hdrg_path = str(pathlib.Path(__file__).resolve().parent.parent / 'herpesdrg-db/herpesdrg-db.tsv')
    hdrg = pd.read_csv(hdrg_path, sep='\t')
    
    co_records = []
    
    # HerpesDRG co-mutations
    hdrg_co = hdrg[hdrg['co_gene'].notna() | hdrg['co_aa'].notna()].copy()
    drugs_list = ['Aciclovir', 'Foscarnet', 'Penciclovir', 'Brivudine', 'Cidofovir', 'Amenamevir', 'Pritelivir']
    
    for idx, r in hdrg_co.iterrows():
        v = r['virus']
        g = r['gene']
        aa = r['aa_change']
        co_g = r['co_gene']
        co_aa = r['co_aa']
        doi = r['ref_doi']
        for d in drugs_list:
            if pd.notna(r[d]):
                co_records.append({
                    'source_database': 'HerpesDRG-DB',
                    'virus': v,
                    'primary_gene': g,
                    'primary_mutation': aa,
                    'co_gene': str(co_g),
                    'co_mutation': str(co_aa),
                    'drug': d,
                    'raw_phenotype': str(r[d]),
                    'reference': str(doi),
                    'exclusion_rationale': 'Multiple genomic alterations present simultaneously; phenotypic resistance cannot be deterministically assigned to a single amino acid site.'
                })
                
    # CMI combos (from df where mutation_category == 'co_occurrence_combo')
    cmi_co = df[df['mutation_category'] == 'co_occurrence_combo']
    for idx, r in cmi_co.iterrows():
        co_records.append({
            'source_database': 'CMI_2025_Review',
            'virus': r['virus'],
            'primary_gene': r['gene'],
            'primary_mutation': r['aa_change'],
            'co_gene': r['gene'],
            'co_mutation': 'Multi-site / Slash notation',
            'drug': r['drug'],
            'raw_phenotype': r['consensus_phenotype'],
            'reference': str(r['sources']),
            'exclusion_rationale': 'Dual/multiple substitution construct; individual site marginal contributions require combinatorial epistasis deconvolution.'
        })
        
    co_df = pd.DataFrame(co_records)
    co_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/co_occurrence_entries.tsv')
    co_df.to_csv(co_path, sep='\t', index=False)
    print(f"\nSaved {len(co_df)} co-occurrence / combo records to: {co_path}")
    
    # 3. Export primary position folds file
    # Strictly HSV-1, primary scope, pure single substitutions, non-conflicted (binary label 0 or 1)
    primary_mask = (
        (df['virus'] == 'HSV-1') &
        (df['is_primary_scope'] == 1) &
        (df['is_single_site'] == 1) &
        (df['target_binary'].isin([0, 1])) &
        (df['is_conflicted'] == 0)
    )
    p_df = df[primary_mask].copy().reset_index(drop=True)
    
    # Create deterministic 5-fold split grouped by position for each pair
    p_df['cv_fold_id'] = -1
    gkf = GroupKFold(n_splits=5)
    
    for (g, d), grp_indices in p_df.groupby(['gene', 'drug']).groups.items():
        sub = p_df.loc[grp_indices]
        for fold_i, (tr_idx, te_idx) in enumerate(gkf.split(sub, groups=sub['pos'])):
            actual_te_indices = grp_indices[te_idx]
            p_df.loc[actual_te_indices, 'cv_fold_id'] = fold_i
            
            # Hard assertion: zero position overlap
            tr_positions = set(sub.iloc[tr_idx]['pos'])
            te_positions = set(sub.iloc[te_idx]['pos'])
            overlap = tr_positions.intersection(te_positions)
            assert len(overlap) == 0, f"LEAKAGE DETECTED in {g}-{d} Fold {fold_i}: overlap={overlap}"
            
    print("\nPrimary Scope Clean Single-Site Dataset Folds:")
    print(p_df.groupby(['gene', 'drug', 'cv_fold_id']).size().unstack(fill_value=0))
    
    folds_out = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/primary_position_folds.tsv')
    cols_to_save = ['virus', 'gene', 'drug', 'aa_change', 'wt', 'pos', 'mut', 'target_binary', 'temporal_split', 'cv_fold_id', 'sources']
    p_df[cols_to_save].to_csv(folds_out, sep='\t', index=False)
    print(f"\nSaved clean primary position folds to: {folds_out}")
    
    # 4. Hash the fold assignments file
    with open(folds_out, 'rb') as f:
        file_bytes = f.read()
        sha256_hash = hashlib.sha256(file_bytes).hexdigest()
        
    hash_log_path = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/fold_assignments_hash.log')
    with open(hash_log_path, 'w', encoding='utf-8') as f:
        f.write("=== REPRODUCIBLE FOLD ASSIGNMENTS SHA-256 INTEGRITY HASH ===\n")
        f.write(f"File: data/processed/primary_position_folds.tsv\n")
        f.write(f"Total Rows: {len(p_df)}\n")
        f.write(f"Primary Scope: HSV-1 UL23-Aciclovir, UL30-Aciclovir, UL30-Foscarnet\n")
        f.write(f"Inclusion Filter: virus=='HSV-1', is_single_site==1, is_conflicted==0, target_binary in [0, 1]\n")
        f.write(f"CV Scheme: 5-Fold Position-Grouped Cross-Validation (GroupKFold)\n")
        f.write(f"Random Seed: 42 (Deterministic GroupKFold)\n")
        f.write(f"SHA-256: {sha256_hash}\n")
        f.write(f"Leakage Verification: All folds asserted zero position overlap between train and test (assert len(overlap) == 0)\n")
        
    print(f"Computed SHA-256: {sha256_hash}")
    print(f"Recorded hash in: {hash_log_path}")

if __name__ == '__main__':
    main()

