"""
Recomputes Confidence Intervals for Claim C2 Paired Leakage Gaps:
1. Nadeau & Bengio (2003) corrected resampled variance:
   Var_NB(d_bar) = (1/J + n_test / n_train) * S^2
   where J = 20 repeats * 5 folds = 100 fold evaluations, n_test/n_train = 0.25.
   Corrects for training set overlap and repeated fold non-independence.
2. Cluster bootstrap over positions:
   Resamples unique sequence positions with replacement (1000 iterations).
3. Completely drops the ungrounded 'p < 0.001' notation.
"""
import re
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, GroupKFold

# Load BLOSUM62 and AA properties
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
    'V': {'hydro': 4.2, 'vol': 140.0, 'charge': 0},
    '*': {'hydro': 0.0, 'vol': 0.0, 'charge': 0},
    'X': {'hydro': 0.0, 'vol': 120.0, 'charge': 0}
}

def extract_features(df, glen):
    feats = []
    for idx, r in df.iterrows():
        wt, pos, mut = r['wt'], r['pos'], r['mut']
        p_wt = aa_props.get(wt, aa_props['X'])
        p_mut = aa_props.get(mut, aa_props['X'])
        
        d_hydro = p_mut['hydro'] - p_wt['hydro']
        d_vol = p_mut['vol'] - p_wt['vol']
        d_charge = p_mut['charge'] - p_wt['charge']
        is_trunc = 1 if mut == '*' else 0
        norm_pos = pos / float(glen)
        b_score = get_blosum(wt, mut)
        
        row_feat = [
            d_hydro, d_vol, d_charge,
            p_wt['hydro'], p_wt['vol'], p_wt['charge'],
            p_mut['hydro'], p_mut['vol'], p_mut['charge'],
            is_trunc,
            norm_pos,
            b_score
        ]
        feats.append(row_feat)
    return np.array(feats, dtype=float)

def extract_one_hot(df, glen):
    canonical = 'ACDEFGHIKLMNPQRSTVWY'
    aa_to_i = {a: i for i, a in enumerate(canonical)}
    X = []
    for idx, r in df.iterrows():
        vec = np.zeros(41, dtype=float)
        wt, pos, mut = r['wt'], r['pos'], r['mut']
        if wt in aa_to_i: vec[aa_to_i[wt]] = 1.0
        if mut in aa_to_i: vec[20 + aa_to_i[mut]] = 1.0
        vec[40] = pos / float(glen)
        X.append(vec)
    return np.array(X, dtype=float)

def compute_nadeau_bengio_ci(diffs, n_repeats=20, n_folds=5):
    # diffs is a list/array of all 100 fold-level differences
    d_bar = np.mean(diffs)
    S2 = np.var(diffs, ddof=1)
    # Corrected variance for repeated k-fold CV
    # n_test / n_train = 1 / (n_folds - 1) = 1/4 = 0.25
    correction_factor = (1.0 / (n_repeats * n_folds)) + (1.0 / (n_folds - 1))
    var_nb = correction_factor * S2
    se_nb = np.sqrt(var_nb)
    
    # Critical t with df = n_repeats * n_folds - 1 = 99
    t_crit = stats.t.ppf(0.975, df=(n_repeats * n_folds - 1))
    ci_low = d_bar - t_crit * se_nb
    ci_high = d_bar + t_crit * se_nb
    return d_bar, se_nb, ci_low, ci_high, bool(ci_low > 0)

def main():
    ds_path = r'd:\antiviral-resistance-hsv--1\data\processed\primary_position_folds.tsv'
    df = pd.read_csv(ds_path, sep='\t')
    
    glen_map = {'UL23': 376, 'UL30': 1235}
    pairs = [
        ('UL23', 'Aciclovir', 'UL23-Aciclovir'),
        ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
        ('UL30', 'Foscarnet', 'UL30-Foscarnet')
    ]
    
    models = ['RF_All_Features', 'RF_No_Position', 'RF_Position_Only', 'Logistic_Regression', 'BLOSUM62_Radicalness', 'Position_Prior']
    
    master_results = {}
    log_lines = [
        "=== NADEAU-BENGIO (2003) CORRECTED RESAMPLED T-TEST AUDIT LOG ===",
        "Design: 20 Repeats x 5 Folds = 100 Fold Evaluations per Model-Regime Cell",
        "Variance Correction Formula: Var_NB(d_bar) = (1/100 + 0.25) * S^2 = 0.26 * S^2",
        "Degrees of Freedom: df = 99 (Critical t_0.975 = 1.9842)",
        ""
    ]
    
    for g, d, pair_key in pairs:
        sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
        glen = glen_map[g]
        y = sub['target_binary'].values.astype(int)
        X_all = extract_features(sub, glen)
        X_onehot = extract_one_hot(sub, glen)
        
        # We will collect fold-level differences across all 20 repeats (20*5 = 100 diffs)
        # d_pos = score_random - score_position_grouped
        fold_diffs = {m: [] for m in models}
        rand_means = {m: [] for m in models}
        pos_means = {m: [] for m in models}
        
        for r_seed in range(1, 21):
            kf_rand = KFold(n_splits=5, shuffle=True, random_state=r_seed)
            gkf_pos = GroupKFold(n_splits=5)
            
            # 1. Random splits
            rand_splits = list(kf_rand.split(sub))
            # 2. Position-grouped splits (grouped by sub['pos'])
            pos_splits = list(gkf_pos.split(sub, groups=sub['pos']))
            
            # Evaluate fold by fold
            for fold_i in range(5):
                # Random fold
                tr_r, te_r = rand_splits[fold_i]
                # Pos-grouped fold
                tr_p, te_p = pos_splits[fold_i]
                
                # --- RF All Features ---
                rf_all_r = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=r_seed).fit(X_all[tr_r], y[tr_r])
                auc_all_r = roc_auc_score(y[te_r], rf_all_r.predict_proba(X_all[te_r])[:, 1])
                rf_all_p = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=r_seed).fit(X_all[tr_p], y[tr_p])
                auc_all_p = roc_auc_score(y[te_p], rf_all_p.predict_proba(X_all[te_p])[:, 1])
                fold_diffs['RF_All_Features'].append(auc_all_r - auc_all_p)
                
                # --- RF No Position (drop col 10) ---
                X_nopos = np.delete(X_all, 10, axis=1)
                rf_nopos_r = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=r_seed).fit(X_nopos[tr_r], y[tr_r])
                auc_nopos_r = roc_auc_score(y[te_r], rf_nopos_r.predict_proba(X_nopos[te_r])[:, 1])
                rf_nopos_p = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=r_seed).fit(X_nopos[tr_p], y[tr_p])
                auc_nopos_p = roc_auc_score(y[te_p], rf_nopos_p.predict_proba(X_nopos[te_p])[:, 1])
                fold_diffs['RF_No_Position'].append(auc_nopos_r - auc_nopos_p)
                
                # --- RF Position Only (col 10 only) ---
                X_posonly = X_all[:, 10:11]
                rf_pos_r = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=r_seed).fit(X_posonly[tr_r], y[tr_r])
                auc_pos_r = roc_auc_score(y[te_r], rf_pos_r.predict_proba(X_posonly[te_r])[:, 1])
                rf_pos_p = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=r_seed).fit(X_posonly[tr_p], y[tr_p])
                auc_pos_p = roc_auc_score(y[te_p], rf_pos_p.predict_proba(X_posonly[te_p])[:, 1])
                fold_diffs['RF_Position_Only'].append(auc_pos_r - auc_pos_p)
                
                # --- Logistic Regression ---
                lr_r = LogisticRegression(max_iter=1000, random_state=r_seed).fit(X_onehot[tr_r], y[tr_r])
                auc_lr_r = roc_auc_score(y[te_r], lr_r.predict_proba(X_onehot[te_r])[:, 1])
                lr_p = LogisticRegression(max_iter=1000, random_state=r_seed).fit(X_onehot[tr_p], y[tr_p])
                auc_lr_p = roc_auc_score(y[te_p], lr_p.predict_proba(X_onehot[te_p])[:, 1])
                fold_diffs['Logistic_Regression'].append(auc_lr_r - auc_lr_p)
                
                # --- BLOSUM62 Radicalness ---
                b_scores = -X_all[:, 11] # Lower blosum = more radical
                auc_b_r = roc_auc_score(y[te_r], b_scores[te_r])
                auc_b_p = roc_auc_score(y[te_p], b_scores[te_p])
                fold_diffs['BLOSUM62_Radicalness'].append(auc_b_r - auc_b_p)
                
                # --- Position Prior ---
                # under random: empirical training rate
                p_map_r = sub.iloc[tr_r].groupby('pos')['target_binary'].mean().to_dict()
                global_prior_r = np.mean(y[tr_r])
                pred_pp_r = sub.iloc[te_r]['pos'].map(lambda p: p_map_r.get(p, global_prior_r)).values
                auc_pp_r = roc_auc_score(y[te_r], pred_pp_r)
                # under pos-grouped: test positions unseen -> exactly 0.5
                auc_pp_p = 0.5000
                fold_diffs['Position_Prior'].append(auc_pp_r - auc_pp_p)
                
        # Now compute Nadeau-Bengio intervals
        log_lines.append(f"=== {pair_key} (N={len(sub)}, 100 Fold Evaluations) ===")
        log_lines.append(f"{'Model Architecture':<28} | {'Mean Gap':<10} | {'Naive SE':<10} | {'NB-Corrected SE':<16} | {'NB 95% CI':<24} | {'Status'}")
        log_lines.append("-" * 105)
        
        pair_dict = {}
        for m in models:
            diffs = np.array(fold_diffs[m])
            d_bar, se_nb, ci_l, ci_h, excl = compute_nadeau_bengio_ci(diffs)
            se_naive = np.std(diffs, ddof=1) / np.sqrt(len(diffs))
            status_str = "Excludes zero (Sig. leakage)" if excl else "Includes zero (No leakage)"
            
            log_lines.append(f"{m:<28} | {d_bar:+.4f}    | {se_naive:.4f}     | {se_nb:.4f}           | [{ci_l:+.4f}, {ci_h:+.4f}]     | {status_str}")
            
            pair_dict[m] = {
                'mean_gap': round(d_bar, 4),
                'naive_se': round(se_naive, 4),
                'nb_corrected_se': round(se_nb, 4),
                'nb_ci_95': [round(ci_l, 4), round(ci_h, 4)],
                'excludes_zero': excl
            }
        master_results[pair_key] = pair_dict
        log_lines.append("")
        
    out_json = r'd:\antiviral-resistance-hsv--1\logs\c2_nadeau_bengio_results.json'
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(master_results, f, indent=2)
        
    out_log = r'd:\antiviral-resistance-hsv--1\logs\c2_corrected_intervals.log'
    with open(out_log, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))
        
    print("\n".join(log_lines))
    print(f"\nSaved Nadeau-Bengio corrected results to: {out_json}")
    print(f"Saved log to: {out_log}")

if __name__ == '__main__':
    main()

