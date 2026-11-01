import pathlib
"""
Temporal Split Confounding & Domain Shift Analysis.
Compares Pre-2016 vs Post-2016 partitions:
1. Class balance (Resistant vs Sensitive)
2. Assay mix
3. Source mix (HerpesDRG vs CMI 2025)
4. Domain classifier: Predicts Pre-2016 vs Post-2016 from NON-LABEL features only.
Reports 5-fold cross-validated AUROC to measure temporal covariate shift.
"""
import re
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

def main():
    ds_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/unified_modeling_dataset.tsv')
    df = pd.read_csv(ds_path, sep='\t')
    
    # Filter to valid temporal split
    df_temp = df[df['temporal_split'].isin(['Pre-2016 (Train)', 'Post-2016 (Test)'])].copy()
    df_temp['is_post_2016'] = (df_temp['temporal_split'] == 'Post-2016 (Test)').astype(int)
    
    log_lines = [
        "=== TEMPORAL SPLIT CONFOUNDING & DOMAIN SHIFT ANALYSIS LOG ===",
        f"Total Samples Analyzed: {len(df_temp)}",
        f"Pre-2016 Samples (Train): {(df_temp['is_post_2016'] == 0).sum()} ({(df_temp['is_post_2016'] == 0).mean()*100:.2f}%)",
        f"Post-2016 Samples (Test): {(df_temp['is_post_2016'] == 1).sum()} ({(df_temp['is_post_2016'] == 1).mean()*100:.2f}%)",
        ""
    ]
    
    # 1. Class Balance Comparison
    pre = df_temp[df_temp['is_post_2016'] == 0]
    post = df_temp[df_temp['is_post_2016'] == 1]
    
    pre_res = (pre['target_binary'] == 1).sum()
    pre_sen = (pre['target_binary'] == 0).sum()
    pre_res_rate = pre_res / (pre_res + pre_sen) if (pre_res + pre_sen) > 0 else 0
    
    post_res = (post['target_binary'] == 1).sum()
    post_sen = (post['target_binary'] == 0).sum()
    post_res_rate = post_res / (post_res + post_sen) if (post_res + post_sen) > 0 else 0
    
    log_lines.append("--- 1. Class Balance Comparison ---")
    log_lines.append(f"Pre-2016:  Resistant = {pre_res}, Sensitive = {pre_sen} -> Resistance Rate = {pre_res_rate*100:.2f}%")
    log_lines.append(f"Post-2016: Resistant = {post_res}, Sensitive = {post_sen} -> Resistance Rate = {post_res_rate*100:.2f}%")
    log_lines.append(f"Class Balance Shift: {post_res_rate - pre_res_rate:+.4f} ({(post_res_rate - pre_res_rate)*100:+.2f} percentage points)")
    log_lines.append("")
    
    # 2. Gene Breakdown Comparison
    log_lines.append("--- 2. Target Gene Breakdown ---")
    for g in df_temp['gene'].unique():
        n_pre = (pre['gene'] == g).sum()
        n_post = (post['gene'] == g).sum()
        log_lines.append(f"Gene {g:<8}: Pre-2016 = {n_pre:<5} ({n_pre/len(pre)*100:5.2f}%) | Post-2016 = {n_post:<5} ({n_post/len(post)*100:5.2f}%)")
    log_lines.append("")
    
    # 3. Drug Breakdown Comparison
    log_lines.append("--- 3. Drug Breakdown ---")
    for d in df_temp['drug'].unique():
        n_pre = (pre['drug'] == d).sum()
        n_post = (post['drug'] == d).sum()
        log_lines.append(f"Drug {d:<12}: Pre-2016 = {n_pre:<5} ({n_pre/len(pre)*100:5.2f}%) | Post-2016 = {n_post:<5} ({n_post/len(post)*100:5.2f}%)")
    log_lines.append("")
    
    # 4. Source Mix (HerpesDRG vs CMI 2025)
    # Check sources string
    pre_cmi = pre['sources'].str.contains('34228537|39772124|40349973|CMI', na=False).sum()
    post_cmi = post['sources'].str.contains('34228537|39772124|40349973|CMI', na=False).sum()
    log_lines.append("--- 4. Literature Source Mix ---")
    log_lines.append(f"Pre-2016 Rows with CMI Review Citations: {pre_cmi} ({pre_cmi/len(pre)*100:.2f}%)")
    log_lines.append(f"Post-2016 Rows with CMI Review Citations: {post_cmi} ({post_cmi/len(post)*100:.2f}%)")
    log_lines.append("")
    
    # 5. Train Domain Classifier (Predicting Pre vs Post from Non-Label Features)
    # Features: gene (one-hot), drug (one-hot), virus (one-hot), normalized position
    glen_map = {'UL23': 376, 'UL30': 1235, 'UL5': 880, 'UL52': 1050, 'UL29': 1200}
    df_temp['norm_pos'] = df_temp.apply(lambda r: r['pos'] / float(glen_map.get(r['gene'], 1000)), axis=1)
    
    X_domain = pd.get_dummies(df_temp[['gene', 'drug', 'virus']], drop_first=True)
    X_domain['norm_pos'] = df_temp['norm_pos'].values
    y_domain = df_temp['is_post_2016'].values
    
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    lr_scores = []
    rf_scores = []
    
    for tr_idx, te_idx in kf.split(X_domain, y_domain):
        X_tr, X_te = X_domain.iloc[tr_idx], X_domain.iloc[te_idx]
        y_tr, y_te = y_domain[tr_idx], y_domain[te_idx]
        
        # Logistic Regression
        lr = LogisticRegression(max_iter=1000, random_state=42)
        lr.fit(X_tr, y_tr)
        p_lr = lr.predict_proba(X_te)[:, 1]
        lr_scores.append(roc_auc_score(y_te, p_lr))
        
        # Random Forest
        rf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
        rf.fit(X_tr, y_tr)
        p_rf = rf.predict_proba(X_te)[:, 1]
        rf_scores.append(roc_auc_score(y_te, p_rf))
        
    log_lines.append("--- 5. Domain Classifier: Predicting Pre-2016 vs Post-2016 from Non-Label Features ---")
    log_lines.append(f"Logistic Regression 5-fold AUROC: {np.mean(lr_scores):.4f} +/- {np.std(lr_scores):.4f} (Folds: {[round(x,3) for x in lr_scores]})")
    log_lines.append(f"Random Forest 5-fold AUROC:       {np.mean(rf_scores):.4f} +/- {np.std(rf_scores):.4f} (Folds: {[round(x,3) for x in rf_scores]})")
    log_lines.append("")
    log_lines.append("Methodological Interpretation:")
    if np.mean(rf_scores) > 0.65:
        log_lines.append(f"  [MODERATE-TO-HIGH DOMAIN SHIFT]: Domain classifier AUROC is {np.mean(rf_scores):.4f} (significantly above 0.50).")
        log_lines.append("  This indicates that post-2016 data exhibits systematic non-random distribution shifts in target drugs, genes, or sequence positions relative to pre-2016 data.")
        log_lines.append("  Therefore, temporal validation performance reflects both temporal generalization AND domain shift across experimental eras.")
    else:
        log_lines.append(f"  [MINIMAL DOMAIN SHIFT]: Domain classifier AUROC is {np.mean(rf_scores):.4f} (near chance).")
        log_lines.append("  Non-label feature distributions are balanced across pre-2016 and post-2016 eras.")
        
    out_log = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/temporal_confounding.log')
    with open(out_log, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))
        
    print("\n".join(log_lines))
    print(f"\nSaved temporal confounding log to: {out_log}")

if __name__ == '__main__':
    main()

