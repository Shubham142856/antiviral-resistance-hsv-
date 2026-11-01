"""
run_spec_exact_tolerance_models.py
Fits the Where x What tolerance-window model using the spec-exact pre-registered Where formula:
  Where(p) = Z(Constraint_PLM) - Z(dist_site)
where Constraint_PLM = - mu_plm_19 (since raw esm2 cache is log P(mut) - log P(wt) <= 0).
Fits all 3 arms with position-clustered bootstrap 95% CIs.
Saves results to logs/where_what_spec_exact_results.json.
"""
import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss, log_loss

# 1. Load data
df = pd.read_csv('data/processed/where_what_features.tsv', sep='\t')

# 2. Compute spec-exact pre-registered features
for g in ['UL23', 'UL30']:
    m = (df['gene'] == g)
    # Z(Constraint): raw mu_plm_19 is negative at conserved positions, so negate it
    z_mu = - (df.loc[m, 'mu_plm_19'] - df.loc[m, 'mu_plm_19'].mean()) / df.loc[m, 'mu_plm_19'].std()
    z_dist = (df.loc[m, 'dist_to_site_angstrom'] - df.loc[m, 'dist_to_site_angstrom'].mean()) / df.loc[m, 'dist_to_site_angstrom'].std()
    df.loc[m, 'where_specexact'] = (z_mu - z_dist) / np.std(z_mu - z_dist)
    # Severity: raw severity_plm is negative for disruptive mutations, so negate it
    df.loc[m, 'sev_specexact'] = - (df.loc[m, 'severity_plm'] - df.loc[m, 'severity_plm'].mean()) / df.loc[m, 'severity_plm'].std()

def fit_model(sub_df, n_boot=1000, seed=42):
    rng = np.random.RandomState(seed)
    W = sub_df['where_specexact'].values
    S = sub_df['sev_specexact'].values
    y = sub_df['target_binary'].values.astype(int)
    
    W_std = (W - np.mean(W)) / np.std(W)
    S_std = (S - np.mean(S)) / np.std(S)
    
    X = np.column_stack([
        W_std,
        S_std,
        S_std**2,
        W_std * S_std,
        W_std * (S_std**2)
    ])
    
    clf = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=seed)
    clf.fit(X, y)
    
    beta_point = [float(clf.intercept_[0])] + [float(c) for c in clf.coef_[0]]
    prob_preds = clf.predict_proba(X)[:, 1]
    auc_fit = float(roc_auc_score(y, prob_preds))
    brier_fit = float(brier_score_loss(y, prob_preds))
    logloss_fit = float(log_loss(y, prob_preds))
    
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
                boot_betas.append([float(clf_b.intercept_[0])] + [float(c) for c in clf_b.coef_[0]])
            except:
                pass
    boot_betas = np.array(boot_betas)
    ci_low = np.percentile(boot_betas, 2.5, axis=0)
    ci_high = np.percentile(boot_betas, 97.5, axis=0)
    
    names = ['Intercept (b0)', 'Where (b1)', 'Severity (b2)', 'Severity^2 (b3)', 'Where*Severity (b4)', 'Where*Severity^2 (b5)']
    coef_dict = {}
    for i, name in enumerate(names):
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

results = {}
targets = [
    ('UL30', 'Aciclovir', 'Primary Hypothesis Target'),
    ('UL30', 'Foscarnet', 'Secondary Target (Power Caveat)'),
    ('UL23', 'Aciclovir', 'Comparative Baseline (Dispensable Target)')
]

print("="*85)
print("SPEC-EXACT 'WHERE x WHAT' TOLERANCE WINDOW TEST RESULTS")
print("="*85)

for gene, drug, role in targets:
    sub = df[(df['gene'] == gene) & (df['drug'] == drug)].copy()
    pair_key = f"{gene}-{drug}"
    pos_cnt = int(sub['target_binary'].sum())
    res = fit_model(sub)
    results[pair_key] = {'role': role, 'N': len(sub), 'positives': pos_cnt, 'spec_exact': res}
    
    print(f"\nTarget: {pair_key} ({role}) | N={len(sub)}, Positives={pos_cnt}")
    print(f"Full Model Fit AUROC: {res['metrics']['auc']:.4f} | Brier: {res['metrics']['brier']:.4f}")
    print("-" * 75)
    for param, pdata in res['coefficients'].items():
        print(f"  {param:<26}: {pdata['point']:+8.4f} | 95% CI: [{pdata['ci_95'][0]:+8.4f}, {pdata['ci_95'][1]:+8.4f}]")

with open('logs/where_what_spec_exact_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nSaved spec-exact results to logs/where_what_spec_exact_results.json")

