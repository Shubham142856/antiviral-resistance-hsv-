import pathlib
import pandas as pd
import numpy as np
import re

def compute_strict_purity():
    db_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/raw/herpesdrg-db/herpesdrg-db.tsv')
    df = pd.read_csv(db_path, sep='\t')
    df_a = df[df['status'] == 'A'].copy()
    df_hsv1 = df_a[df_a['virus'] == 'HSV1'].copy()
    
    def extract_pos(aa):
        m = re.search(r'\d+', str(aa))
        return int(m.group()) if m else np.nan
    df_hsv1['pos'] = df_hsv1['aa_change'].apply(extract_pos)
    
    def classify(val, th):
        if pd.isna(val): return np.nan
        s = str(val).strip()
        if s.lower() == 'polymorphism': return 'Polymorphism'
        if s.lower() == 'resistant': return 'Resistant'
        try:
            fc = float(s)
            return 'Resistant' if fc >= th else 'Susceptible'
        except:
            return 'Unknown'
            
    pairs = [
        ('UL23', 'Aciclovir', 'Sauerbrei 2016 / Castor 2012 (FC >= 2.0 / 3.0)'),
        ('UL30', 'Aciclovir', 'Sauerbrei 2016 / HerpesDRG 2024 (FC >= 2.0 / 3.0)'),
        ('UL30', 'Foscarnet', 'Sauerbrei 2016 / Burrel 2013 (FC >= 2.0 / 3.0)')
    ]
    
    print("==================================================================")
    print("STRICT LABEL PURITY AUDIT (Restricted to Multi-Substitution Positions)")
    print("==================================================================")
    
    for gene, drug, cite in pairs:
        sub = df_hsv1[(df_hsv1['gene'] == gene) & (df_hsv1[drug].notna())].copy()
        sub['label_2x'] = sub[drug].apply(lambda x: classify(x, 2.0))
        sub['label_3x'] = sub[drug].apply(lambda x: classify(x, 3.0))
        
        # Binary target: 1 = Resistant, 0 = Susceptible/Polymorphism
        sub['target_2x'] = sub['label_2x'].map({'Resistant': 1, 'Susceptible': 0, 'Polymorphism': 0})
        sub['target_3x'] = sub['label_3x'].map({'Resistant': 1, 'Susceptible': 0, 'Polymorphism': 0})
        sub = sub[sub['target_2x'].notna()].copy()
        
        # Group by position
        pos_stats = sub.groupby('pos').agg(
            n_rows=('aa_change', 'count'),
            n_unique_muts=('aa_change', 'nunique'),
            n_res_2x=('target_2x', 'sum'),
            mean_res_2x=('target_2x', 'mean'),
            n_res_3x=('target_3x', 'sum'),
            mean_res_3x=('target_3x', 'mean')
        ).reset_index()
        
        total_pos = len(pos_stats)
        single_sub_pos = pos_stats[pos_stats['n_unique_muts'] == 1]
        multi_sub_pos = pos_stats[pos_stats['n_unique_muts'] >= 2]
        
        # Multi-substitution purity at 2.0x
        multi_homo_2x = multi_sub_pos[(multi_sub_pos['mean_res_2x'] == 0.0) | (multi_sub_pos['mean_res_2x'] == 1.0)]
        multi_mixed_2x = multi_sub_pos[(multi_sub_pos['mean_res_2x'] > 0.0) & (multi_sub_pos['mean_res_2x'] < 1.0)]
        
        # Multi-substitution purity at 3.0x
        multi_homo_3x = multi_sub_pos[(multi_sub_pos['mean_res_3x'] == 0.0) | (multi_sub_pos['mean_res_3x'] == 1.0)]
        multi_mixed_3x = multi_sub_pos[(multi_sub_pos['mean_res_3x'] > 0.0) & (multi_sub_pos['mean_res_3x'] < 1.0)]
        
        print(f"\n--- {gene} - {drug} (HSV-1) [{cite}] ---")
        print(f"Total evaluated positions: {total_pos}")
        print(f"  Single-substitution positions: {len(single_sub_pos)} ({len(single_sub_pos)/total_pos*100:.1f}%) [Trivially homogeneous]")
        print(f"  Multi-substitution positions (>=2 distinct mutations): {len(multi_sub_pos)} ({len(multi_sub_pos)/total_pos*100:.1f}%)")
        print(f"  Among multi-substitution positions (Threshold >= 2.0x):")
        print(f"    - Homogeneous (all resistant OR all susceptible): {len(multi_homo_2x)} ({len(multi_homo_2x)/len(multi_sub_pos)*100:.1f}%)")
        print(f"    - Heterogeneous / Mixed (different mutations have different labels): {len(multi_mixed_2x)} ({len(multi_mixed_2x)/len(multi_sub_pos)*100:.1f}%)")
        print(f"  Among multi-substitution positions (Threshold >= 3.0x):")
        print(f"    - Homogeneous: {len(multi_homo_3x)} ({len(multi_homo_3x)/len(multi_sub_pos)*100:.1f}%)")
        print(f"    - Heterogeneous / Mixed: {len(multi_mixed_3x)} ({len(multi_mixed_3x)/len(multi_sub_pos)*100:.1f}%)")
        
        if len(multi_mixed_2x) > 0:
            print("  Example heterogeneous positions (different substitutions -> different outcomes):")
            for idx, r in multi_mixed_2x.head(5).iterrows():
                p = int(r['pos'])
                muts = sub[sub['pos'] == p][['aa_change', drug, 'label_2x']].to_dict(orient='records')
                print(f"    Pos {p}: {[(m['aa_change'], m['label_2x']) for m in muts]}")

if __name__ == '__main__':
    compute_strict_purity()
