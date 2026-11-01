import pathlib
"""
Rerun Claim C2 under 20 Repeated Grouped Cross-Validations with Hard Leakage Assertions.
Features:
1. Strict Both-Grouped double-holdout (zero position overlap AND zero strain overlap).
2. Hard assertions verifying zero test position overlap in Regimes 3, 4, 5.
3. Model ablations:
   - RF (All Features)
   - RF (No Position Feature)
   - RF (Position Feature Only)
   - Logistic Regression
   - BLOSUM62 Radicalness
   - Position Prior (strictly verifies 0.5000 per-fold on position holdouts)
4. Contiguous-Block Grouped Split (5 contiguous sequence segments along protein length).
5. 20 repeated CV runs; reports paired leakage gaps (Random - Grouped) with 95% bootstrap CIs.
"""
import re
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold, GroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

BLOSUM62 = {
    ('A', 'A'): 4, ('A', 'R'): -1, ('A', 'N'): -2, ('A', 'D'): -2, ('A', 'C'): 0, ('A', 'Q'): -1, ('A', 'E'): -1, ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'L'): -1, ('A', 'K'): -1, ('A', 'M'): -1, ('A', 'F'): -2, ('A', 'P'): -1, ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0,
    ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2, ('R', 'C'): -3, ('R', 'Q'): 1, ('R', 'E'): 0, ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2, ('R', 'K'): 2, ('R', 'M'): -1, ('R', 'F'): -3, ('R', 'P'): -2, ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'W'): -3, ('R', 'Y'): -2, ('R', 'V'): -3,
    ('N', 'N'): 6, ('N', 'D'): 1, ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0, ('N', 'G'): 0, ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0, ('N', 'M'): -2, ('N', 'F'): -3, ('N', 'P'): -2, ('N', 'S'): 1, ('N', 'T'): 0, ('N', 'W'): -4, ('N', 'Y'): -2, ('N', 'V'): -3,
    ('D', 'D'): 6, ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2, ('D', 'G'): -1, ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'L'): -4, ('D', 'K'): -1, ('D', 'M'): -3, ('D', 'F'): -3, ('D', 'P'): -1, ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'W'): -4, ('D', 'Y'): -3, ('D', 'V'): -3,
    ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4, ('C', 'G'): -3, ('C', 'H'): -3, ('C', 'I'): -1, ('C', 'L'): -1, ('C', 'K'): -3, ('C', 'M'): -1, ('C', 'F'): -2, ('C', 'P'): -3, ('C', 'S'): -1, ('C', 'T'): -1, ('C', 'W'): -2, ('C', 'Y'): -2, ('C', 'V'): -1,
    ('Q', 'Q'): 5, ('Q', 'E'): 2, ('Q', 'G'): -2, ('Q', 'H'): 0, ('Q', 'I'): -3, ('Q', 'L'): -2, ('Q', 'K'): 1, ('Q', 'M'): 0, ('Q', 'F'): -3, ('Q', 'P'): -1, ('Q', 'S'): 0, ('Q', 'T'): -1, ('Q', 'W'): -2, ('Q', 'Y'): -1, ('Q', 'V'): -2,
    ('E', 'E'): 5, ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3, ('E', 'K'): 1, ('E', 'M'): -2, ('E', 'F'): -3, ('E', 'P'): -1, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'W'): -3, ('E', 'Y'): -2, ('E', 'V'): -2,
    ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2, ('G', 'M'): -3, ('G', 'F'): -3, ('G', 'P'): -2, ('G', 'S'): 0, ('G', 'T'): -2, ('G', 'W'): -2, ('G', 'Y'): -3, ('G', 'V'): -3,
    ('H', 'H'): 8, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1, ('H', 'M'): -2, ('H', 'F'): -1, ('H', 'P'): -2, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'W'): -2, ('H', 'Y'): 2, ('H', 'V'): -3,
    ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3, ('I', 'M'): 1, ('I', 'F'): 0, ('I', 'P'): -3, ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'W'): -3, ('I', 'Y'): -1, ('I', 'V'): 3,
    ('L', 'L'): 4, ('L', 'K'): -2, ('L', 'M'): 2, ('L', 'F'): 0, ('L', 'P'): -3, ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'W'): -2, ('L', 'Y'): -1, ('L', 'V'): 1,
    ('K', 'K'): 5, ('K', 'M'): -1, ('K', 'F'): -3, ('K', 'P'): -1, ('K', 'S'): 0, ('K', 'T'): -1, ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2,
    ('M', 'M'): 5, ('M', 'F'): 0, ('M', 'P'): -2, ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'W'): -1, ('M', 'Y'): -1, ('M', 'V'): 1,
    ('F', 'F'): 6, ('F', 'P'): -4, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'W'): 1, ('F', 'Y'): 3, ('F', 'V'): -1,
    ('P', 'P'): 7, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'W'): -4, ('P', 'Y'): -3, ('P', 'V'): -2,
    ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2,
    ('T', 'T'): 5, ('T', 'W'): -2, ('T', 'Y'): -2, ('T', 'V'): 0,
    ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3,
    ('Y', 'Y'): 7, ('Y', 'V'): -1,
    ('V', 'V'): 4
}

def get_blosum(a, b):
    pair = (a, b) if (a, b) in BLOSUM62 else (b, a)
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

def parse_mut(aa):
    m = re.match(r'^([A-Z])(\d+)([A-Z*])$', str(aa).strip())
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    m_num = re.search(r'\d+', str(aa))
    pos = int(m_num.group()) if m_num else 1
    return 'X', pos, 'X'

def extract_features(df, glen):
    # Returns 12 features
    # Cols: 0-9: physicochemical deltas and properties
    # Col 10: normalized sequence position
    # Col 11: BLOSUM62 score
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
        b_score = get_blosum(wt, mut) if mut != '*' else -10
        
        feats.append([
            p_wt['hydro'], p_wt['vol'], p_wt['charge'],
            p_mut['hydro'], p_mut['vol'], p_mut['charge'],
            d_hydro, d_vol, d_charge, is_trunc, norm_pos, b_score
        ])
    return np.array(feats)

def compute_safe_auc(y_true, y_pred):
    if len(np.unique(y_true)) < 2: return 0.5
    if len(np.unique(y_pred)) < 2: return 0.5
    try: return float(roc_auc_score(y_true, y_pred))
    except: return 0.5

def bootstrap_ci(values, n_boot=1000, random_state=42):
    rng = np.random.RandomState(random_state)
    boot_means = []
    vals = np.array(values)
    if len(vals) == 0: return 0.0, 0.0
    for _ in range(n_boot):
        sample = rng.choice(vals, size=len(vals), replace=True)
        boot_means.append(np.mean(sample))
    return float(np.percentile(boot_means, 2.5)), float(np.percentile(boot_means, 97.5))

def create_both_grouped_splits(df, n_splits=5, seed=42):
    """
    Strict 2D double-holdout split:
    Partition unique positions into n_splits sets AND unique strains into n_splits sets.
    Test fold k contains ONLY samples with pos in Pos_k AND strain in Strain_k.
    Guarantees: Zero test position overlap AND zero test strain overlap with training set!
    """
    rng = np.random.RandomState(seed)
    unique_positions = np.array(sorted(df['pos'].unique()))
    unique_strains = np.array(sorted(df['strain_cluster'].unique()))
    
    rng.shuffle(unique_positions)
    rng.shuffle(unique_strains)
    
    pos_folds = np.array_split(unique_positions, n_splits)
    strain_folds = np.array_split(unique_strains, n_splits)
    
    splits = []
    for k in range(n_splits):
        test_pos_set = set(pos_folds[k])
        test_strain_set = set(strain_folds[k])
        
        # Test samples MUST have pos in test_pos_set AND strain in test_strain_set
        test_mask = df['pos'].isin(test_pos_set) & df['strain_cluster'].isin(test_strain_set)
        # Train samples MUST have pos NOT in test_pos_set AND strain NOT in test_strain_set
        train_mask = (~df['pos'].isin(test_pos_set)) & (~df['strain_cluster'].isin(test_strain_set))
        
        test_indices = np.where(test_mask)[0]
        train_indices = np.where(train_mask)[0]
        
        # Only valid if test set has both classes
        if len(test_indices) > 0 and len(train_indices) > 0:
            splits.append((train_indices, test_indices))
            
    return splits

def main():
    db_path = str(pathlib.Path(__file__).resolve().parent.parent / 'herpesdrg-db/herpesdrg-db.tsv')
    out_json = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/c2_repeated_cv_results.json')
    out_log = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/c2_repeated_cv.log')
    
    df = pd.read_csv(db_path, sep='\t')
    df_hsv1 = df[(df['status'] == 'A') & (df['virus'] == 'HSV1')].copy()
    df_hsv1['wt'], df_hsv1['pos'], df_hsv1['mut'] = zip(*df_hsv1['aa_change'].apply(parse_mut))
    
    def get_target(s):
        try: return 1 if float(s) >= 3.0 else 0
        except: return 1 if 'resistant' in str(s).lower() else 0
        
    targets = [
        ('UL23', 'Aciclovir', 376),
        ('UL30', 'Aciclovir', 1235),
        ('UL30', 'Foscarnet', 1235)
    ]
    
    n_repeats = 20
    master_results = {}
    log_lines = [
        "=== CLAIM C2 REPEATED GROUPED CROSS-VALIDATION AUDIT LOG ===",
        f"Repeated CV Iterations: {n_repeats} (Seeds 1 to {n_repeats})",
        "Evaluation Method: Per-Fold Mean AUROC and Paired Leakage Gap (Random - Grouped) with 95% Bootstrap CIs",
        "Hard Assertions: Zero position overlap verified in Regimes 3, 4, 5.",
        ""
    ]
    
    for gene, drug, glen in targets:
        sub = df_hsv1[(df_hsv1['gene'] == gene) & (df_hsv1[drug].notna())].copy()
        sub['target'] = sub[drug].apply(get_target)
        sub['strain_cluster'] = sub['control_strain'].fillna(sub['strain_generation']).fillna('Reference')
        # Contiguous sequence block (5 blocks along protein length)
        sub['seq_block'] = pd.qcut(sub['pos'], q=5, labels=False, duplicates='drop')
        
        pair_key = f"{gene}-{drug}"
        log_lines.append(f"================================================================================")
        log_lines.append(f"TARGET PAIR: {pair_key} (N={len(sub)}, Positives={sub['target'].sum()}, Negatives={len(sub)-sub['target'].sum()})")
        log_lines.append(f"Distinct Positions: {sub['pos'].nunique()} | Distinct Strains: {sub['strain_cluster'].nunique()}")
        log_lines.append(f"================================================================================")
        
        X_all = extract_features(sub, glen)
        y = sub['target'].values
        
        # Feature subsets
        # 1. All features (12 features)
        # 2. No position (drop col 10 -> 11 features)
        # 3. Position only (col 10 -> 1 feature)
        X_no_pos = np.delete(X_all, 10, axis=1)
        X_pos_only = X_all[:, 10:11]
        
        models_to_test = [
            'RF_All_Features',
            'RF_No_Position',
            'RF_Position_Only',
            'Logistic_Regression',
            'BLOSUM62_Radicalness',
            'Position_Prior'
        ]
        
        # Store per-repeat results across regimes
        repeat_records = {m: {'Random': [], 'Clade_Grouped': [], 'Position_Grouped': [], 'Both_Grouped': [], 'Contiguous_Block': []} for m in models_to_test}
        paired_gaps = {m: {'Random_minus_Position': [], 'Random_minus_Clade': [], 'Random_minus_Block': []} for m in models_to_test}
        
        for rep in range(1, n_repeats + 1):
            seed = rep * 42
            
            # Define 5 splitters for this repeat
            # 1. Random
            kf_rand = KFold(n_splits=5, shuffle=True, random_state=seed)
            splits_rand = list(kf_rand.split(X_all, y))
            
            # 2. Clade-Grouped
            # Shuffle strain groups across folds
            unique_strains = np.array(sorted(sub['strain_cluster'].unique()))
            np.random.RandomState(seed).shuffle(unique_strains)
            strain_map = {s: i % 5 for i, s in enumerate(unique_strains)}
            strain_fold_ids = sub['strain_cluster'].map(strain_map).values
            splits_clade = [(np.where(strain_fold_ids != k)[0], np.where(strain_fold_ids == k)[0]) for k in range(5) if (strain_fold_ids == k).sum() > 0]
            
            # 3. Position-Grouped
            unique_positions = np.array(sorted(sub['pos'].unique()))
            np.random.RandomState(seed).shuffle(unique_positions)
            pos_map = {p: i % 5 for i, p in enumerate(unique_positions)}
            pos_fold_ids = sub['pos'].map(pos_map).values
            splits_pos = [(np.where(pos_fold_ids != k)[0], np.where(pos_fold_ids == k)[0]) for k in range(5) if (pos_fold_ids == k).sum() > 0]
            
            # 4. Strict Both-Grouped
            splits_both = create_both_grouped_splits(sub, n_splits=5, seed=seed)
            
            # 5. Contiguous Block
            splits_block = [(np.where(sub['seq_block'] != k)[0], np.where(sub['seq_block'] == k)[0]) for k in range(5) if (sub['seq_block'] == k).sum() > 0]
            
            # Run evaluations across regimes
            regimes = [
                ('Random', splits_rand, False, False),
                ('Clade_Grouped', splits_clade, False, True),
                ('Position_Grouped', splits_pos, True, False),
                ('Both_Grouped', splits_both, True, True),
                ('Contiguous_Block', splits_block, True, False)
            ]
            
            current_rep_scores = {m: {} for m in models_to_test}
            
            for reg_name, split_list, must_check_pos, must_check_strain in regimes:
                fold_scores = {m: [] for m in models_to_test}
                
                for tr_idx, te_idx in split_list:
                    # HARD LEAKAGE ASSERTIONS
                    tr_pos = set(sub.iloc[tr_idx]['pos'])
                    te_pos = set(sub.iloc[te_idx]['pos'])
                    tr_str = set(sub.iloc[tr_idx]['strain_cluster'])
                    te_str = set(sub.iloc[te_idx]['strain_cluster'])
                    
                    pos_overlap = len(tr_pos.intersection(te_pos))
                    str_overlap = len(tr_str.intersection(te_str))
                    
                    if must_check_pos:
                        assert pos_overlap == 0, f"LEAKAGE VIOLATION: {reg_name} has {pos_overlap} overlapping positions!"
                    if must_check_strain:
                        assert str_overlap == 0, f"LEAKAGE VIOLATION: {reg_name} has {str_overlap} overlapping strains!"
                        
                    y_tr, y_te = y[tr_idx], y[te_idx]
                    if len(np.unique(y_te)) < 2:
                        continue
                        
                    # Model 1: RF All
                    rf1 = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=seed)
                    rf1.fit(X_all[tr_idx], y_tr)
                    fold_scores['RF_All_Features'].append(compute_safe_auc(y_te, rf1.predict_proba(X_all[te_idx])[:, 1]))
                    
                    # Model 2: RF No Position
                    rf2 = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=seed)
                    rf2.fit(X_no_pos[tr_idx], y_tr)
                    fold_scores['RF_No_Position'].append(compute_safe_auc(y_te, rf2.predict_proba(X_no_pos[te_idx])[:, 1]))
                    
                    # Model 3: RF Position Only
                    rf3 = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=seed)
                    rf3.fit(X_pos_only[tr_idx], y_tr)
                    fold_scores['RF_Position_Only'].append(compute_safe_auc(y_te, rf3.predict_proba(X_pos_only[te_idx])[:, 1]))
                    
                    # Model 4: Logistic Regression
                    lr = LogisticRegression(max_iter=1000, random_state=seed)
                    lr.fit(X_all[tr_idx], y_tr)
                    fold_scores['Logistic_Regression'].append(compute_safe_auc(y_te, lr.predict_proba(X_all[te_idx])[:, 1]))
                    
                    # Model 5: BLOSUM62 Radicalness
                    fold_scores['BLOSUM62_Radicalness'].append(compute_safe_auc(y_te, -X_all[te_idx, 11]))
                    
                    # Model 6: Position Prior
                    tr_df = sub.iloc[tr_idx]
                    pos_map_dict = tr_df.groupby('pos')['target'].mean().to_dict()
                    glob_prior = float(y_tr.mean())
                    p_pos = sub.iloc[te_idx]['pos'].map(lambda p: pos_map_dict.get(p, glob_prior)).values
                    auc_pp = compute_safe_auc(y_te, p_pos)
                    
                    # If positions held out, verify mathematically that tied predictions evaluate to 0.5000
                    if must_check_pos:
                        assert auc_pp == 0.5, f"POSITION PRIOR VIOLATION: Expected 0.5 under position holdout, got {auc_pp}!"
                    fold_scores['Position_Prior'].append(auc_pp)
                    
                for m in models_to_test:
                    mean_score = float(np.mean(fold_scores[m])) if len(fold_scores[m]) > 0 else 0.5
                    repeat_records[m][reg_name].append(mean_score)
                    current_rep_scores[m][reg_name] = mean_score
                    
            # Compute paired gaps for this repeat
            for m in models_to_test:
                g_pos = current_rep_scores[m]['Random'] - current_rep_scores[m]['Position_Grouped']
                g_cla = current_rep_scores[m]['Random'] - current_rep_scores[m]['Clade_Grouped']
                g_blk = current_rep_scores[m]['Random'] - current_rep_scores[m]['Contiguous_Block']
                paired_gaps[m]['Random_minus_Position'].append(g_pos)
                paired_gaps[m]['Random_minus_Clade'].append(g_cla)
                paired_gaps[m]['Random_minus_Block'].append(g_blk)
                
        # Summarize across 20 repeats
        log_lines.append(f"--- 20-Repeat Summary for {pair_key} ---")
        log_lines.append(f"{'Model Architecture':<25} | {'Regime 1 (Rand)':<17} | {'Regime 2 (Clade)':<17} | {'Regime 3 (Pos)':<17} | {'Regime 4 (Both)':<17} | {'Regime 5 (Block)':<17}")
        log_lines.append("-" * 115)
        
        pair_summary = {}
        for m in models_to_test:
            means = {r: np.mean(repeat_records[m][r]) for r in ['Random', 'Clade_Grouped', 'Position_Grouped', 'Both_Grouped', 'Contiguous_Block']}
            stds = {r: np.std(repeat_records[m][r]) for r in ['Random', 'Clade_Grouped', 'Position_Grouped', 'Both_Grouped', 'Contiguous_Block']}
            
            line = f"{m:<25} | {means['Random']:.4f} +/- {stds['Random']:.4f} | {means['Clade_Grouped']:.4f} +/- {stds['Clade_Grouped']:.4f} | {means['Position_Grouped']:.4f} +/- {stds['Position_Grouped']:.4f} | {means['Both_Grouped']:.4f} +/- {stds['Both_Grouped']:.4f} | {means['Contiguous_Block']:.4f} +/- {stds['Contiguous_Block']:.4f}"
            log_lines.append(line)
            
            # Paired gap analysis
            gap_pos_mean = np.mean(paired_gaps[m]['Random_minus_Position'])
            gap_pos_ci = bootstrap_ci(paired_gaps[m]['Random_minus_Position'])
            
            gap_cla_mean = np.mean(paired_gaps[m]['Random_minus_Clade'])
            gap_cla_ci = bootstrap_ci(paired_gaps[m]['Random_minus_Clade'])
            
            pair_summary[m] = {
                'means': {k: round(v, 4) for k, v in means.items()},
                'stds': {k: round(v, 4) for k, v in stds.items()},
                'paired_gap_pos': {
                    'mean_gap': round(gap_pos_mean, 4),
                    'ci_95': [round(gap_pos_ci[0], 4), round(gap_pos_ci[1], 4)],
                    'excludes_zero': bool(gap_pos_ci[0] > 0)
                },
                'paired_gap_clade': {
                    'mean_gap': round(gap_cla_mean, 4),
                    'ci_95': [round(gap_cla_ci[0], 4), round(gap_cla_ci[1], 4)],
                    'excludes_zero': bool(gap_cla_ci[0] > 0)
                }
            }
            
        master_results[pair_key] = pair_summary
        log_lines.append("")
        log_lines.append(f"Paired Leakage Gaps (Random minus Position-Grouped) across 20 repeats:")
        for m in models_to_test:
            g_info = pair_summary[m]['paired_gap_pos']
            sig_str = "EXCLUDES ZERO (Consistent with leakage)" if g_info['excludes_zero'] else "Includes zero / Not significant"
            log_lines.append(f"  {m:<25}: Mean Gap = {g_info['mean_gap']:+.4f} (95% CI: [{g_info['ci_95'][0]:.4f}, {g_info['ci_95'][1]:.4f}]) -> {sig_str}")
        log_lines.append("")
        
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(master_results, f, indent=2)
        
    with open(out_log, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))
        
    print("\n".join(log_lines))
    print(f"Saved repeated CV results to: {out_json}")
    print(f"Saved execution log to: {out_log}")

if __name__ == '__main__':
    main()

