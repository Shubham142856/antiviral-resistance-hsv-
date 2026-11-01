import pathlib
import pandas as pd
import numpy as np
import re
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold, GroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

# Standard physicochemical scales (Kyte-Doolittle hydropathy, molecular volume, net charge)
aa_properties = {
    'A': {'hydro': 1.8, 'vol': 88.6, 'charge': 0},
    'R': {'hydro': -4.5, 'vol': 173.4, 'charge': 1},
    'N': {'hydro': -3.5, 'vol': 114.1, 'charge': 0},
    'D': {'hydro': -3.5, 'vol': 111.1, 'charge': -1},
    'C': {'hydro': 2.5, 'vol': 108.5, 'charge': 0},
    'Q': {'hydro': -3.5, 'vol': 143.8, 'charge': 0},
    'E': {'hydro': -3.5, 'vol': 138.4, 'charge': -1},
    'G': {'hydro': -0.4, 'vol': 60.1, 'charge': 0},
    'H': {'hydro': -3.2, 'vol': 153.2, 'charge': 0},
    'I': {'hydro': 4.5, 'vol': 166.7, 'charge': 0},
    'L': {'hydro': 3.8, 'vol': 166.7, 'charge': 0},
    'K': {'hydro': -3.9, 'vol': 168.6, 'charge': 1},
    'M': {'hydro': 1.9, 'vol': 162.9, 'charge': 0},
    'F': {'hydro': 2.8, 'vol': 189.9, 'charge': 0},
    'P': {'hydro': -1.6, 'vol': 112.7, 'charge': 0},
    'S': {'hydro': -0.8, 'vol': 89.0, 'charge': 0},
    'T': {'hydro': -0.7, 'vol': 116.1, 'charge': 0},
    'W': {'hydro': -0.9, 'vol': 227.8, 'charge': 0},
    'Y': {'hydro': -1.3, 'vol': 193.6, 'charge': 0},
    'V': {'hydro': 4.2, 'vol': 140.0, 'charge': 0},
    '*': {'hydro': 0.0, 'vol': 0.0, 'charge': 0},
    'X': {'hydro': 0.0, 'vol': 120.0, 'charge': 0}
}

blosum62_diag = {'A':4,'R':5,'N':6,'D':6,'C':9,'Q':5,'E':5,'G':6,'H':8,'I':4,'L':4,'K':5,'M':5,'F':6,'P':7,'S':4,'T':5,'W':11,'Y':7,'V':4}

def parse_mut(aa):
    m = re.match(r'^([A-Z])(\d+)([A-Z*])$', str(aa).strip())
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    m_num = re.search(r'\d+', str(aa))
    pos = int(m_num.group()) if m_num else 1
    return 'X', pos, 'X'

def extract_features(df):
    feats = []
    for idx, r in df.iterrows():
        wt, pos, mut = r['wt'], r['pos'], r['mut']
        p_wt = aa_properties.get(wt, aa_properties['X'])
        p_mut = aa_properties.get(mut, aa_properties['X'])
        
        # Features:
        # 1. Delta properties
        d_hydro = p_mut['hydro'] - p_wt['hydro']
        d_vol = p_mut['vol'] - p_wt['vol']
        d_charge = p_mut['charge'] - p_wt['charge']
        is_trunc = 1 if mut == '*' else 0
        
        # 2. Position coordinate normalized
        norm_pos = pos / 376.0 # for TK
        
        # 3. Known functional domains for TK (ATP binding: 50-65; Nucleoside binding: 160-180)
        in_atp = 1 if (50 <= pos <= 65) else 0
        in_nucleoside = 1 if (160 <= pos <= 180) else 0
        
        feats.append([
            p_wt['hydro'], p_wt['vol'], p_wt['charge'],
            p_mut['hydro'], p_mut['vol'], p_mut['charge'],
            d_hydro, d_vol, d_charge, is_trunc,
            norm_pos, in_atp, in_nucleoside
        ])
def extract_features_gene(df, gene_len):
    feats = []
    for idx, r in df.iterrows():
        wt, pos, mut = r['wt'], r['pos'], r['mut']
        p_wt = aa_properties.get(wt, aa_properties['X'])
        p_mut = aa_properties.get(mut, aa_properties['X'])
        
        d_hydro = p_mut['hydro'] - p_wt['hydro']
        d_vol = p_mut['vol'] - p_wt['vol']
        d_charge = p_mut['charge'] - p_wt['charge']
        is_trunc = 1 if mut == '*' else 0
        norm_pos = pos / float(gene_len)
        
        feats.append([
            p_wt['hydro'], p_wt['vol'], p_wt['charge'],
            p_mut['hydro'], p_mut['vol'], p_mut['charge'],
            d_hydro, d_vol, d_charge, is_trunc, norm_pos
        ])
    return np.array(feats)

def bootstrap_auc_ci(y_true, y_pred, n_boot=1000, random_state=42):
    rng = np.random.RandomState(random_state)
    boot_aucs = []
    indices = np.arange(len(y_true))
    for _ in range(n_boot):
        sample_idx = rng.choice(indices, size=len(indices), replace=True)
        if len(np.unique(y_true[sample_idx])) < 2:
            continue
        try:
            boot_aucs.append(roc_auc_score(y_true[sample_idx], y_pred[sample_idx]))
        except:
            continue
    if len(boot_aucs) == 0:
        return 0.5, 0.5
    return np.percentile(boot_aucs, 2.5), np.percentile(boot_aucs, 97.5)

def run_c2_test():
    import json
    db_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/raw/herpesdrg-db/herpesdrg-db.tsv')
    df = pd.read_csv(db_path, sep='\t')
    df_a = df[df['status'] == 'A'].copy()
    df_hsv1 = df_a[df_a['virus'] == 'HSV1'].copy()
    df_hsv1['wt'], df_hsv1['pos'], df_hsv1['mut'] = zip(*df_hsv1['aa_change'].apply(parse_mut))
    
    targets = [
        ('UL23', 'Aciclovir', 376),
        ('UL30', 'Aciclovir', 1235),
        ('UL30', 'Foscarnet', 1235)
    ]
    
    def classify(x):
        s = str(x).strip().lower()
        if s == 'polymorphism': return 0
        if s == 'resistant': return 1
        try: return 1 if float(s) >= 2.0 else 0
        except: return np.nan
        
    all_results = {}
    
    for gene, drug, glen in targets:
        sub = df_hsv1[(df_hsv1['gene'] == gene) & (df_hsv1[drug].notna())].copy()
        sub['target'] = sub[drug].apply(classify)
        sub = sub[sub['target'].notna()].copy()
        sub['target'] = sub['target'].astype(int)
        
        # Strain cluster
        sub['strain_cluster'] = sub['control_strain'].fillna(sub['strain_generation']).fillna('Reference')
        sub['pos_cluster'] = sub['pos'].astype(str)
        sub['pos_strain_cluster'] = sub['pos'].astype(str) + "_" + sub['strain_cluster']
        
        X = extract_features_gene(sub, glen)
        y = sub['target'].values
        
        pair_key = f"{gene}-{drug}"
        print(f"\n=== Evaluating Claim C2 on {pair_key} (N={len(sub)}, Pos={y.sum()}, Neg={len(y)-y.sum()}) ===")
        
        models = {
            'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42),
            'Random Forest': RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
        }
        
        regimes = [
            ('Regime 1: Random Split', KFold(n_splits=5, shuffle=True, random_state=42), None),
            ('Regime 2: Position-Grouped Split', GroupKFold(n_splits=5), sub['pos']),
            ('Regime 3: Strain-Grouped Split', GroupKFold(n_splits=5), sub['strain_cluster']),
            ('Regime 4: Position+Strain Grouped', GroupKFold(n_splits=5), sub['pos_strain_cluster'])
        ]
        
        pair_res = []
        for m_name, model in models.items():
            for r_name, splitter, groups in regimes:
                preds = np.zeros(len(y))
                if groups is None:
                    splits = list(splitter.split(X, y))
                else:
                    splits = list(splitter.split(X, y, groups=groups))
                    
                for train_idx, test_idx in splits:
                    X_train, y_train = X[train_idx], y[train_idx]
                    X_test, y_test = X[test_idx], y[test_idx]
                    
                    if len(np.unique(y_train)) < 2:
                        preds[test_idx] = float(np.mean(y_train))
                        continue
                    model.fit(X_train, y_train)
                    preds[test_idx] = model.predict_proba(X_test)[:, 1]
                    
                auc = roc_auc_score(y, preds)
                ci_low, ci_high = bootstrap_auc_ci(y, preds)
                pair_res.append({
                    'pair': pair_key,
                    'model': m_name,
                    'regime': r_name,
                    'auroc': round(float(auc), 4),
                    'ci_95': [round(float(ci_low), 4), round(float(ci_high), 4)]
                })
                print(f"[{m_name}] {r_name}: AUROC = {auc:.4f} (95% CI: [{ci_low:.4f}, {ci_high:.4f}])")
                
        all_results[pair_key] = pair_res
        
    out_file = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/c2_split_evaluation.json')
    with open(out_file, 'w', encoding='utf-8') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved detailed evaluation to {out_file}")

if __name__ == '__main__':
    run_c2_test()
