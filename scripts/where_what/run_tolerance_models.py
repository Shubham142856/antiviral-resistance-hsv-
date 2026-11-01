"""
run_tolerance_models.py
Fits the pre-registered Where x What tolerance-window model:
  logit(P(Resistant)) = b0 + b1*Where + b2*Sev + b3*Sev^2 + b4*(Where*Sev) + b5*(Where*Sev^2)
Computes position-clustered bootstrap 95% CIs and evaluates explicit falsifiers.
Anchored to pre-registration: docs/WHERE_WHAT_PREREGISTRATION.md
"""

import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss, log_loss

# 1. Load data
df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')
print(f"Loaded {len(df)} variants across {df['gene'].nunique()} targets.")

targets = [
    ('UL30', 'Aciclovir', 'Primary Hypothesis Target'),
    ('UL30', 'Foscarnet', 'Secondary Target (Power Caveat)'),
    ('UL23', 'Aciclovir', 'Comparative Baseline (Dispensable Target)')
]

def fit_tolerance_model(sub_df, sev_col, where_col='where_composite', n_boot=1000, seed=42):
    rng = np.random.RandomState(seed)
    
    # Feature construction
    W = sub_df[where_col].values
    S = sub_df[sev_col].values
    y = sub_df['target_binary'].values
    
    # Standardize predictors for numerical stability and interpretable coefficients
    W_std = (W - np.mean(W)) / np.std(W)
    S_std = (S - np.mean(S)) / np.std(S)
    
    # Design matrix: [W, S, S^2, W*S, W*S^2]
    X = np.column_stack([
        W_std,
        S_std,
        S_std**2,
        W_std * S_std,
        W_std * (S_std**2)
    ])
    
    # Fit point estimate with L2 penalization (C=1.0)
    clf = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=seed)
    clf.fit(X, y)
    
    beta_point = [float(clf.intercept_[0])] + [float(c) for c in clf.coef_[0]]
    
    # Performance metrics
    prob_preds = clf.predict_proba(X)[:, 1]
    auc_fit = float(roc_auc_score(y, prob_preds))
    brier_fit = float(brier_score_loss(y, prob_preds))
    logloss_fit = float(log_loss(y, prob_preds))
    
    # Position-clustered bootstrap for 95% CIs
    unique_positions = np.array(sorted(sub_df['pos'].unique()))
    pos_indices = {p: sub_df[sub_df['pos'] == p].index.values for p in unique_positions}
    
    boot_betas = []
    for _ in range(n_boot):
        sample_pos = rng.choice(unique_positions, size=len(unique_positions), replace=True)
        sample_idx = np.concatenate([pos_indices[p] for p in sample_pos])
        
        y_b = y[sub_df.index.get_indexer_for(sample_idx)]
        X_b = X[sub_df.index.get_indexer_for(sample_idx)]
        
        if len(np.unique(y_b)) > 1:
            try:
                clf_b = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=500)
                clf_b.fit(X_b, y_b)
                b_b = [float(clf_b.intercept_[0])] + [float(c) for c in clf_b.coef_[0]]
                boot_betas.append(b_b)
            except Exception:
                continue
                
    boot_betas = np.array(boot_betas)
    ci_low = np.percentile(boot_betas, 2.5, axis=0)
    ci_high = np.percentile(boot_betas, 97.5, axis=0)
    
    feature_names = ['Intercept (b0)', 'Where (b1)', 'Severity (b2)', 'Severity^2 (b3)', 'Where*Severity (b4)', 'Where*Severity^2 (b5)']
    
    coef_dict = {}
    for i, name in enumerate(feature_names):
        coef_dict[name] = {
            'point': beta_point[i],
            'ci_95': [float(ci_low[i]), float(ci_high[i])]
        }
        
    return {
        'coefficients': coef_dict,
        'metrics': {
            'auc': auc_fit,
            'brier': brier_fit,
            'log_loss': logloss_fit
        }
    }

results_master = {}

print("\n" + "="*85)
print("PRE-REGISTERED 'WHERE x WHAT' TOLERANCE WINDOW TEST RESULTS")
print("="*85)

for gene, drug, role in targets:
    sub = df[(df['gene'] == gene) & (df['drug'] == drug)].copy()
    pair_key = f"{gene}-{drug}"
    pos_cnt = int(sub['target_binary'].sum())
    
    print(f"\nTarget: {pair_key} ({role}) | N={len(sub)}, Positives={pos_cnt} ({(pos_cnt/len(sub))*100:.1f}%)")
    print("-" * 85)
    
    results_master[pair_key] = {'role': role, 'N': len(sub), 'positives': pos_cnt, 'models': {}}
    
    for sev_name, sev_col in [('PLM Relative Severity', 'z_severity_plm'), ('BLOSUM62 Radicalness', 'z_severity_blosum')]:
        res = fit_tolerance_model(sub, sev_col)
        results_master[pair_key]['models'][sev_name] = res
        
        print(f"\n--- Severity Metric: {sev_name} ---")
        print(f"{'Parameter':<28} | {'Point Est':<10} | {'95% Cluster Boot CI'}")
        print("-" * 65)
        for param, pdata in res['coefficients'].items():
            print(f"{param:<28} | {pdata['point']:<10.4f} | [{pdata['ci_95'][0]:.4f}, {pdata['ci_95'][1]:.4f}]")
            
        print(f"Model Fit: AUROC = {res['metrics']['auc']:.4f} | Brier = {res['metrics']['brier']:.4f} | LogLoss = {res['metrics']['log_loss']:.4f}")

# 2. Evaluate Pre-registered Falsification Rules
print("\n" + "="*85)
print("EVALUATION OF PRE-REGISTERED FALSIFICATION CRITERIA")
print("="*85)

# A. UL30-ACV Tolerance Window Falsifier: b5 < 0 (negative interaction curvature)
ul30_acv_res = results_master['UL30-Aciclovir']['models']['PLM Relative Severity']
b5_point = ul30_acv_res['coefficients']['Where*Severity^2 (b5)']['point']
b5_ci = ul30_acv_res['coefficients']['Where*Severity^2 (b5)']['ci_95']

print(f"1. Primary Target (UL30-Aciclovir) Interaction Curvature (b5):")
print(f"   Point Estimate: {b5_point:.4f} | 95% CI: [{b5_ci[0]:.4f}, {b5_ci[1]:.4f}]")
if b5_ci[1] < 0:
    print("   -> OUTCOME: PREDICTION CONFIRMED! 95% CI strictly excludes zero in negative direction (b5 < 0).")
    print("      Confirms tolerance window: at constrained sites (Where > 0), intermediate severity peaks resistance.")
elif b5_ci[0] <= 0 <= b5_ci[1]:
    if b5_point < 0:
        print("   -> OUTCOME: DIRECTIONAL SUPPORT ONLY. Point estimate is negative (b5 < 0), but 95% CI includes zero.")
        print("      Tolerance window hypothesis is NOT statistically definitive under strict falsifier.")
    else:
        print("   -> OUTCOME: FALSIFIED. Point estimate is non-negative and 95% CI spans zero.")
else:
    print("   -> OUTCOME: FALSIFIED. 95% CI is strictly positive.")

# B. UL23-ACV Monotonicity Falsifier: b3 ~ 0, b5 ~ 0, b2 > 0
ul23_acv_res = results_master['UL23-Aciclovir']['models']['PLM Relative Severity']
b2_point_ul23 = ul23_acv_res['coefficients']['Severity (b2)']['point']
b2_ci_ul23 = ul23_acv_res['coefficients']['Severity (b2)']['ci_95']
b3_point_ul23 = ul23_acv_res['coefficients']['Severity^2 (b3)']['point']
b3_ci_ul23 = ul23_acv_res['coefficients']['Severity^2 (b3)']['ci_95']
b5_point_ul23 = ul23_acv_res['coefficients']['Where*Severity^2 (b5)']['point']
b5_ci_ul23 = ul23_acv_res['coefficients']['Where*Severity^2 (b5)']['ci_95']

print(f"\n2. Dispensable Activator (UL23-Aciclovir) Monotonicity Test:")
print(f"   Severity Linear (b2): {b2_point_ul23:.4f} | 95% CI: [{b2_ci_ul23[0]:.4f}, {b2_ci_ul23[1]:.4f}]")
print(f"   Severity Curvature (b3): {b3_point_ul23:.4f} | 95% CI: [{b3_ci_ul23[0]:.4f}, {b3_ci_ul23[1]:.4f}]")
print(f"   Where*Sev^2 (b5): {b5_point_ul23:.4f} | 95% CI: [{b5_ci_ul23[0]:.4f}, {b5_ci_ul23[1]:.4f}]")

if b2_ci_ul23[0] > 0 and b3_ci_ul23[0] <= 0 <= b3_ci_ul23[1]:
    print("   -> OUTCOME: PREDICTION CONFIRMED! Linear severity is strictly positive, curvature is non-significant.")
    print("      Confirms monotonic loss-of-function on UL23.")
else:
    print("   -> OUTCOME: Review coefficient intervals above.")

# Save master results
with open('logs/where_what_test_results.json', 'w') as f:
    json.dump(results_master, f, indent=2)

print("\n" + "="*85)
print("Complete results saved to logs/where_what_test_results.json")
print("="*85)
