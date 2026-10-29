import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold

BLOSUM62 = {
    ('A','A'): 4, ('A','R'): -1, ('A','N'): -2, ('A','D'): -2, ('A','C'): 0, ('A','Q'): -1, ('A','E'): -1, ('A','G'): 0, ('A','H'): -2, ('A','I'): -1, ('A','L'): -1, ('A','K'): -1, ('A','M'): -1, ('A','F'): -2, ('A','P'): -1, ('A','S'): 1, ('A','T'): 0, ('A','W'): -3, ('A','Y'): -2, ('A','V'): 0,
    ('R','R'): 5, ('R','N'): 0, ('R','D'): -2, ('R','C'): -3, ('R','Q'): 1, ('R','E'): 0, ('R','G'): -2, ('R','H'): 0, ('R','I'): -3, ('R','L'): -2, ('R','K'): 2, ('R','M'): -1, ('R','F'): -3, ('R','P'): -2, ('R','S'): -1, ('R','T'): -1, ('R','W'): -3, ('R','Y'): -2, ('R','V'): -3,
    ('N','N'): 6, ('N','D'): 1, ('N','C'): -3, ('N','Q'): 0, ('N','E'): 0, ('N','G'): 0, ('N','H'): 1, ('N','I'): -3, ('N','L'): -3, ('N','K'): 0, ('N','M'): -2, ('N','F'): -3, ('N','P'): -2, ('N','S'): 1, ('N','T'): 0, ('N','W'): -4, ('N','Y'): -2, ('N','V'): -3,
    ('D','D'): 6, ('D','C'): -3, ('D','Q'): 0, ('D','E'): 2, ('D','G'): -1, ('D','H'): -1, ('D','I'): -3, ('D','L'): -4, ('D','K'): -1, ('D','M'): -3, ('D','F'): -3, ('D','P'): -1, ('D','S'): 0, ('D','T'): -1, ('D','W'): -4, ('D','Y'): -3, ('D','V'): -3,
    ('C','C'): 9, ('C','Q'): -3, ('C','E'): -4, ('C','G'): -3, ('C','H'): -3, ('C','I'): -1, ('C','L'): -1, ('C','K'): -3, ('C','M'): -1, ('C','F'): -2, ('C','P'): -3, ('C','S'): -1, ('C','T'): -1, ('C','W'): -2, ('C','Y'): -2, ('C','V'): -1,
    ('Q','Q'): 5, ('Q','E'): 2, ('Q','G'): -2, ('Q','H'): 0, ('Q','I'): -3, ('Q','L'): -2, ('Q','K'): 1, ('Q','M'): 0, ('Q','F'): -3, ('Q','P'): -1, ('Q','S'): 0, ('Q','T'): -1, ('Q','W'): -2, ('Q','Y'): -1, ('Q','V'): -2,
    ('E','E'): 5, ('E','G'): -2, ('E','H'): 0, ('E','I'): -3, ('E','L'): -3, ('E','K'): 1, ('E','M'): -2, ('E','F'): -3, ('E','P'): -1, ('E','S'): 0, ('E','T'): -1, ('E','W'): -3, ('E','Y'): -2, ('E','V'): -2,
    ('G','G'): 6, ('G','H'): -2, ('G','I'): -4, ('G','L'): -4, ('G','K'): -2, ('G','M'): -3, ('G','F'): -3, ('G','P'): -2, ('G','S'): 0, ('G','T'): -2, ('G','W'): -2, ('G','Y'): -3, ('G','V'): -3,
    ('H','H'): 8, ('H','I'): -3, ('H','L'): -3, ('H','K'): -1, ('H','M'): -2, ('H','F'): -1, ('H','P'): -2, ('H','S'): -1, ('H','T'): -2, ('H','W'): -2, ('H','Y'): 2, ('H','V'): -3,
    ('I','I'): 4, ('I','L'): 2, ('I','K'): -3, ('I','M'): 1, ('I','F'): 0, ('I','P'): -3, ('I','S'): -2, ('I','T'): -1, ('I','W'): -3, ('I','Y'): -1, ('I','V'): 3,
    ('L','L'): 4, ('L','K'): -2, ('L','M'): 2, ('L','F'): 0, ('L','P'): -3, ('L','S'): -2, ('L','T'): -1, ('L','W'): -2, ('L','Y'): -1, ('L','V'): 1,
    ('K','K'): 5, ('K','M'): -1, ('K','F'): -3, ('K','P'): -1, ('K','S'): 0, ('K','T'): -1, ('K','W'): -3, ('K','Y'): -2, ('K','V'): -2,
    ('M','M'): 5, ('M','F'): 0, ('M','P'): -2, ('M','S'): -1, ('M','T'): -1, ('M','W'): -1, ('M','Y'): -1, ('M','V'): 1,
    ('F','F'): 6, ('F','P'): -4, ('F','S'): -2, ('F','T'): -2, ('F','W'): 1, ('F','Y'): 3, ('F','V'): -1,
    ('P','P'): 7, ('P','S'): -1, ('P','T'): -1, ('P','W'): -4, ('P','Y'): -3, ('P','V'): -2,
    ('S','S'): 4, ('S','T'): 1, ('S','W'): -3, ('S','Y'): -2, ('S','V'): -2,
    ('T','T'): 5, ('T','W'): -2, ('T','Y'): -2, ('T','V'): 0,
    ('W','W'): 11, ('W','Y'): 2, ('W','V'): -3,
    ('Y','Y'): 7, ('Y','V'): -1,
    ('V','V'): 4
}

def get_blosum(a1, a2):
    pair = (a1, a2) if (a1, a2) in BLOSUM62 else (a2, a1)
    return BLOSUM62.get(pair, -4)

aa_props = {
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
    'V': {'hydro': 4.2, 'vol': 140.0, 'charge': 0}
}

def extract_features(df, glen):
    feats = []
    for idx, r in df.iterrows():
        wt, mut, pos = r['wt'], r['mut'], r['pos']
        p_wt = aa_props.get(wt, {'hydro': 0, 'vol': 120, 'charge': 0})
        p_mut = aa_props.get(mut, {'hydro': 0, 'vol': 120, 'charge': 0})
        b_score = get_blosum(wt, mut)
        row_feat = [
            p_wt['hydro'], p_wt['vol'], p_wt['charge'],
            p_mut['hydro'], p_mut['vol'], p_mut['charge'],
            p_mut['hydro'] - p_wt['hydro'],
            p_mut['vol'] - p_wt['vol'],
            p_mut['charge'] - p_wt['charge'],
            abs(p_mut['charge'] - p_wt['charge']),
            pos / glen,
            b_score
        ]
        feats.append(row_feat)
    return np.array(feats)

def run_cluster_bootstrap():
    ds_path = r'd:\antiviral-resistance-hsv--1\data\processed\primary_position_folds.tsv'
    df = pd.read_csv(ds_path, sep='\t')
    glen_map = {'UL23': 376.0, 'UL30': 1235.0}
    
    pairs = [
        ('UL23', 'Aciclovir', 'UL23-Aciclovir'),
        ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
        ('UL30', 'Foscarnet', 'UL30-Foscarnet')
    ]
    
    models = ['RF_All_Features', 'RF_No_Position', 'RF_Position_Only', 'Logistic_Regression', 'BLOSUM62_Radicalness', 'Position_Prior']
    n_boot = 1000
    rng = np.random.RandomState(42)
    
    results = {}
    log_lines = [
        "=== CLAIM C2 CLUSTER BOOTSTRAP OVER POSITIONS AUDIT LOG ===",
        f"Bootstrap Iterations: {n_boot} (Resampling Unique Sequence Positions with Replacement)",
        "Metric: Paired Leakage Gap (Naive Random AUROC - Position-Grouped AUROC)",
        ""
    ]
    
    for g, d, pair_key in pairs:
        sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
        glen = glen_map[g]
        y = sub['target_binary'].values.astype(int)
        X_all = extract_features(sub, glen)
        positions = np.array(sorted(sub['pos'].unique()))
        n_pos = len(positions)
        
        # Precompute out-of-fold predictions on original dataset for:
        # 1. Random K-fold (5 folds)
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        preds_rand = {m: np.zeros(len(sub)) for m in models}
        for tr, te in kf.split(sub):
            # RF All
            rf_all = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42).fit(X_all[tr], y[tr])
            preds_rand['RF_All_Features'][te] = rf_all.predict_proba(X_all[te])[:, 1]
            # RF No Pos
            X_nopos = np.delete(X_all, 10, axis=1)
            rf_nopos = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42).fit(X_nopos[tr], y[tr])
            preds_rand['RF_No_Position'][te] = rf_nopos.predict_proba(X_nopos[te])[:, 1]
            # RF Pos Only
            X_posonly = X_all[:, 10:11]
            rf_pos = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42).fit(X_posonly[tr], y[tr])
            preds_rand['RF_Position_Only'][te] = rf_pos.predict_proba(X_posonly[te])[:, 1]
            # LR
            lr = LogisticRegression(max_iter=1000, random_state=42).fit(X_all[tr], y[tr])
            preds_rand['Logistic_Regression'][te] = lr.predict_proba(X_all[te])[:, 1]
            # BLOSUM
            preds_rand['BLOSUM62_Radicalness'][te] = -X_all[te, 11]
            # Position Prior
            p_map = sub.iloc[tr].groupby('pos')['target_binary'].mean().to_dict()
            g_prior = np.mean(y[tr])
            preds_rand['Position_Prior'][te] = sub.iloc[te]['pos'].map(lambda p: p_map.get(p, g_prior)).values
            
        # 2. Position-grouped K-fold (using cv_fold_id in primary_position_folds.tsv)
        preds_pos = {m: np.zeros(len(sub)) for m in models}
        for fold_i in range(5):
            tr = np.where(sub['cv_fold_id'] != fold_i)[0]
            te = np.where(sub['cv_fold_id'] == fold_i)[0]
            # RF All
            rf_all = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42).fit(X_all[tr], y[tr])
            preds_pos['RF_All_Features'][te] = rf_all.predict_proba(X_all[te])[:, 1]
            # RF No Pos
            X_nopos = np.delete(X_all, 10, axis=1)
            rf_nopos = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42).fit(X_nopos[tr], y[tr])
            preds_pos['RF_No_Position'][te] = rf_nopos.predict_proba(X_nopos[te])[:, 1]
            # RF Pos Only
            X_posonly = X_all[:, 10:11]
            rf_pos = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42).fit(X_posonly[tr], y[tr])
            preds_pos['RF_Position_Only'][te] = rf_pos.predict_proba(X_posonly[te])[:, 1]
            # LR
            lr = LogisticRegression(max_iter=1000, random_state=42).fit(X_all[tr], y[tr])
            preds_pos['Logistic_Regression'][te] = lr.predict_proba(X_all[te])[:, 1]
            # BLOSUM
            preds_pos['BLOSUM62_Radicalness'][te] = -X_all[te, 11]
            # Position Prior
            preds_pos['Position_Prior'][te] = 0.5000
            
        # Cluster bootstrap over unique positions
        pos_to_indices = sub.groupby('pos').indices
        boot_diffs = {m: [] for m in models}
        
        for b in range(n_boot):
            sampled_positions = rng.choice(positions, size=n_pos, replace=True)
            sampled_idx = []
            for p in sampled_positions:
                sampled_idx.extend(pos_to_indices[p])
            sampled_idx = np.array(sampled_idx)
            
            y_b = y[sampled_idx]
            if len(np.unique(y_b)) < 2:
                continue
                
            for m in models:
                auc_r = roc_auc_score(y_b, preds_rand[m][sampled_idx])
                # For position prior under pos-grouped, tied predictions -> 0.5
                if m == 'Position_Prior':
                    auc_p = 0.5000
                else:
                    auc_p = roc_auc_score(y_b, preds_pos[m][sampled_idx])
                boot_diffs[m].append(auc_r - auc_p)
                
        log_lines.append(f"=== {pair_key} (N={len(sub)}, {n_pos} Unique Positions, {len(boot_diffs['RF_All_Features'])} Valid Bootstraps) ===")
        log_lines.append(f"{'Model Architecture':<28} | {'Mean Gap':<10} | {'Bootstrap SE':<12} | {'95% Position-Cluster CI':<26} | {'Verdict'}")
        log_lines.append("-" * 105)
        
        pair_dict = {}
        for m in models:
            arr = np.array(boot_diffs[m])
            mean_gap = float(np.mean(arr))
            se_boot = float(np.std(arr, ddof=1))
            ci_low = float(np.percentile(arr, 2.5))
            ci_high = float(np.percentile(arr, 97.5))
            excl_zero = bool(ci_low > 0)
            v_str = "EXCLUDES ZERO (Significant leakage)" if excl_zero else "Spans zero (No leakage)"
            
            log_lines.append(f"{m:<28} | {mean_gap:+.4f}    | {se_boot:.4f}       | [{ci_low:+.4f}, {ci_high:+.4f}]         | {v_str}")
            pair_dict[m] = {
                'mean_gap': round(mean_gap, 4),
                'se_boot': round(se_boot, 4),
                'ci_95_cluster_bootstrap': [round(ci_low, 4), round(ci_high, 4)],
                'excludes_zero': excl_zero
            }
        results[pair_key] = pair_dict
        log_lines.append("")
        
    out_json = r'd:\antiviral-resistance-hsv--1\logs\c2_position_cluster_bootstrap_results.json'
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
        
    out_log = r'd:\antiviral-resistance-hsv--1\logs\c2_position_cluster_bootstrap.log'
    with open(out_log, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))
        
    print("\n".join(log_lines))

if __name__ == '__main__':
    run_cluster_bootstrap()
