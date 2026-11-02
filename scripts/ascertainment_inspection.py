import pandas as pd
import re

herpes = pd.read_csv('herpesdrg-db/herpesdrg-db.tsv', sep='\t')
hsv1 = herpes[herpes['virus'] == 'HSV1'].copy()

def parse_mut(s):
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', str(s).strip())
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    return None, None, None

parsed = [parse_mut(x) for x in hsv1['aa_change']]
hsv1['wt'] = [p[0] for p in parsed]
hsv1['pos'] = [p[1] for p in parsed]
hsv1['mut'] = [p[2] for p in parsed]
hsv1['gene_std'] = hsv1['gene'].replace({'TK': 'UL23', 'POL': 'UL30', 'DPOL': 'UL30'})

feat = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')

for g, d in [('UL30', 'Aciclovir'), ('UL30', 'Foscarnet'), ('UL23', 'Aciclovir')]:
    sub_feat = feat[(feat['gene'] == g) & (feat['drug'] == d)]
    m = pd.merge(sub_feat, hsv1, left_on=['gene', 'pos', 'wt', 'mut'], right_on=['gene_std', 'pos', 'wt', 'mut'], how='inner')
    
    def is_num(val):
        try:
            float(str(val).strip())
            return True
        except:
            return False
            
    m['has_numeric_fc'] = [is_num(x) for x in m[d]]
    num_df = m[m['has_numeric_fc']]
    n_res = int(num_df['target_binary'].sum())
    n_sens = int(len(num_df) - n_res)
    print(f"=== {g}-{d} ===")
    print(f"Total merged: {len(m)}")
    print(f"Has numeric measured fold-change: N={len(num_df)} (Resistant={n_res}, Sensitive={n_sens})")
    print("Top strain_generation categories:")
    for sg, count in m['strain_generation'].value_counts().head(5).items():
        sub_sg = m[m['strain_generation'] == sg]
        print(f"  {sg}: N={len(sub_sg)}, Pos={int(sub_sg['target_binary'].sum())}")
