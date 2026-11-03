"""
run_split_degradation_benchmark.py
==================================
Confirmatory execution of the Feature-Family x Split-Scheme Degradation Benchmark.

Anchored strictly to protocol:
  docs/SPLIT_DEGRADATION_CONFIRMATORY_PROTOCOL.md (Version 1.0)

Evaluates:
  - Three Targets: UL30-Aciclovir (N=332), UL30-Foscarnet (N=311), UL23-Aciclovir (N=286)
  - Three Split Schemes:
      1. Position-Grouped 5-Fold (`cv_fold_id`)
      2. 3D k-Means Spatial 5-Fold (`spatial_cluster_id`)
      3. Contiguous Sequence Block 5-Fold (`seq_block_id`)
  - Four Supervised Feature Families:
      1. Family A: Positional Identity (`norm_pos`)
      2. Family B: Evolutionary Conservation (`blosum62`, `neg_entropy`, `weighted_cons`, `pssm`)
      3. Family C: Protein Language Models (`esm2_650m`, `protbert`, `esm2_35m`)
      4. Family D: Multi-Modal Combined (All features)
  - Explicit Negative Control:
      Fixed Zero-Shot BLOSUM62 (unfitted, invariant by construction)
  - Models:
      1. RandomForestClassifier (n_estimators=100, max_depth=5, random_state=42)
      2. LogisticRegression (penalty='l2', C=1.0, random_state=42)
  - Resampling:
      B = 2,000 Paired Bootstrap Resamples (seed 42)

Output:
  - logs/split_degradation_benchmark_results.json
"""

import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
N_BOOT = 2000
RNG_SEED = 42

PROTEIN_LENGTHS = {
    'UL30': 1235,
    'UL23': 376
}

FEATURE_FAMILIES = {
    'Positional Identity': ['norm_pos'],
    'Evolutionary Conservation': ['score_blosum62', 'score_neg_entropy', 'score_weighted_cons', 'score_pssm'],
    'Protein Language Models': ['score_esm2_650m', 'score_protbert', 'score_esm2_35m'],
    'Multi-Modal Combined': ['norm_pos', 'score_blosum62', 'score_neg_entropy', 'score_weighted_cons', 'score_pssm', 'score_esm2_650m', 'score_protbert', 'score_esm2_35m']
}

SPLIT_SCHEMES = [
    ('Position-Grouped', 'cv_fold_id'),
    ('3D k-Means Spatial', 'spatial_cluster_id'),
    ('Contiguous Block', 'seq_block_id')
]

TARGETS = [
    ('UL30-Aciclovir', 'UL30', 'Aciclovir'),
    ('UL30-Foscarnet', 'UL30', 'Foscarnet'),
    ('UL23-Aciclovir', 'UL23', 'Aciclovir')
]

def run_oof_cv(df_sub, feat_cols, split_col, model_type='rf'):
    y = df_sub['target_binary'].values
    oof_preds = np.zeros(len(df_sub))
    folds = sorted(df_sub[split_col].unique())
    
    X_raw = df_sub[feat_cols].values
    
    for f in folds:
        train_idx = (df_sub[split_col] != f).values
        test_idx = (df_sub[split_col] == f).values
        
        X_train, y_train = X_raw[train_idx], y[train_idx]
        X_test = X_raw[test_idx]
        
        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        if model_type == 'rf':
            clf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
            clf.fit(X_train_scaled, y_train)
            oof_preds[test_idx] = clf.predict_proba(X_test_scaled)[:, 1]
        elif model_type == 'lr':
            clf = LogisticRegression(penalty='l2', C=1.0, random_state=42, max_iter=1000)
            clf.fit(X_train_scaled, y_train)
            oof_preds[test_idx] = clf.predict_proba(X_test_scaled)[:, 1]
            
    auc = float(roc_auc_score(y, oof_preds))
    return auc, oof_preds

def main():
    print("=" * 85)
    print("FEATURE-FAMILY x SPLIT-SCHEME DEGRADATION BENCHMARK")
    print("Protocol: docs/SPLIT_DEGRADATION_CONFIRMATORY_PROTOCOL.md")
    print("=" * 85)

    data_path = os.path.join(ROOT, 'data', 'processed', 'primary_3d_spatial_folds.tsv')
    df = pd.read_csv(data_path, sep='\t')

    # Add normalized position
    df['norm_pos'] = [r['pos'] / PROTEIN_LENGTHS[r['gene']] for _, r in df.iterrows()]

    results = {
        'protocol': {
            'document': 'docs/SPLIT_DEGRADATION_CONFIRMATORY_PROTOCOL.md',
            'version': '1.0',
            'status': 'CONFIRMATORY EXECUTION',
            'n_boot': N_BOOT,
            'seed': RNG_SEED
        },
        'targets': {}
    }

    rng = np.random.RandomState(RNG_SEED)

    for target_name, gene, drug in TARGETS:
        print(f"\nEvaluating Target: {target_name} ({gene} - {drug})")
        sub = df[(df['gene'] == gene) & (df['drug'] == drug)].copy().reset_index(drop=True)
        
        # Subset with 3D coordinates for spatial holdout comparison
        sub_3d = sub[sub['has_3d_coord'] == 1].copy().reset_index(drop=True)
        print(f"  Cohort: N = {len(sub)} (3D subset N = {len(sub_3d)}), {int(sub['target_binary'].sum())} R / {int((sub['target_binary']==0).sum())} S")

        target_record = {
            'n_total': len(sub),
            'n_with_3d': len(sub_3d),
            'n_resistant': int(sub['target_binary'].sum()),
            'n_sensitive': int((sub['target_binary'] == 0).sum()),
            'models': {'random_forest': {}, 'logistic_regression': {}},
            'negative_control_zero_shot_blosum': {}
        }

        # 1. Negative Control: Unfitted Zero-Shot BLOSUM62 across splits
        y_all = sub['target_binary'].values
        base_blosum = sub['score_blosum62'].values
        auc_blosum = float(roc_auc_score(y_all, base_blosum))

        target_record['negative_control_zero_shot_blosum'] = {
            'position_grouped_auc': round(auc_blosum, 4),
            'contiguous_block_auc': round(auc_blosum, 4),
            'degradation_delta': 0.0000,
            'status': 'VERIFIED INVARIANT BY CONSTRUCTION (Delta = 0.0000)'
        }

        # 2. Supervised Models across Feature Families & Split Schemes
        for m_type, m_name in [('rf', 'random_forest'), ('lr', 'logistic_regression')]:
            print(f"\n  Model: {m_name.upper()}")
            print(f"    %-32s %-20s %8s %8s %8s" % ("Feature Family", "Position-Grouped", "3D Spatial", "Contig-Block", "Delta(PG - Block)"))
            print("    " + "-" * 82)

            target_record['models'][m_name] = {}

            for fam_name, feat_cols in FEATURE_FAMILIES.items():
                # Full cohort on Position-Grouped and Contiguous-Block
                auc_pg, preds_pg = run_oof_cv(sub, feat_cols, 'cv_fold_id', model_type=m_type)
                auc_block, preds_block = run_oof_cv(sub, feat_cols, 'seq_block_id', model_type=m_type)

                # 3D spatial subset
                auc_pg_3d, preds_pg_3d = run_oof_cv(sub_3d, feat_cols, 'cv_fold_id', model_type=m_type)
                auc_spatial, preds_spatial = run_oof_cv(sub_3d, feat_cols, 'spatial_cluster_id', model_type=m_type)
                auc_block_3d, preds_block_3d = run_oof_cv(sub_3d, feat_cols, 'seq_block_id', model_type=m_type)

                delta_degrade = auc_pg - auc_block
                delta_spatial = auc_pg_3d - auc_spatial

                # Paired Bootstrap on OOF Predictions
                boot_deltas = []
                boot_auc_pg = []
                boot_auc_block = []

                n_variants = len(sub)
                y_sub = sub['target_binary'].values

                for b in range(N_BOOT):
                    boot_idx = rng.choice(n_variants, size=n_variants, replace=True)
                    y_b = y_sub[boot_idx]
                    if len(np.unique(y_b)) < 2:
                        continue
                    a_pg = roc_auc_score(y_b, preds_pg[boot_idx])
                    a_bl = roc_auc_score(y_b, preds_block[boot_idx])
                    boot_auc_pg.append(a_pg)
                    boot_auc_block.append(a_bl)
                    boot_deltas.append(a_pg - a_bl)

                ci_pg = [float(np.percentile(boot_auc_pg, 2.5)), float(np.percentile(boot_auc_pg, 97.5))]
                ci_block = [float(np.percentile(boot_auc_block, 2.5)), float(np.percentile(boot_auc_block, 97.5))]
                ci_delta = [float(np.percentile(boot_deltas, 2.5)), float(np.percentile(boot_deltas, 97.5))]
                p_val_degrade = float(np.mean(np.array(boot_deltas) <= 0))

                print(f"    %-32s %8.4f [%.3f, %.3f]  %8.4f  %8.4f [%.3f, %.3f]  %+.4f [%+.3f, %+.3f] (p=%.4f)"
                      % (fam_name, auc_pg, ci_pg[0], ci_pg[1], auc_spatial, auc_block, ci_block[0], ci_block[1], delta_degrade, ci_delta[0], ci_delta[1], p_val_degrade))

                target_record['models'][m_name][fam_name] = {
                    'features': feat_cols,
                    'position_grouped_auc': round(auc_pg, 4),
                    'position_grouped_ci_95': [round(ci_pg[0], 4), round(ci_pg[1], 4)],
                    'spatial_3d_auc': round(auc_spatial, 4),
                    'contiguous_block_auc': round(auc_block, 4),
                    'contiguous_block_ci_95': [round(ci_block[0], 4), round(ci_block[1], 4)],
                    'degradation_delta_pg_minus_block': round(delta_degrade, 4),
                    'degradation_delta_ci_95': [round(ci_delta[0], 4), round(ci_delta[1], 4)],
                    'degradation_p_value': round(p_val_degrade, 4),
                    'spatial_subset_comparison': {
                        'n_3d': len(sub_3d),
                        'position_grouped_auc': round(auc_pg_3d, 4),
                        'spatial_3d_auc': round(auc_spatial, 4),
                        'contiguous_block_auc': round(auc_block_3d, 4),
                        'spatial_degradation_delta': round(delta_spatial, 4)
                    }
                }

        results['targets'][target_name] = target_record

    # Save to logs
    out_path = os.path.join(ROOT, 'logs', 'split_degradation_benchmark_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)

    print("\n" + "=" * 85)
    print(f"Wrote complete degradation benchmark results to: {os.path.relpath(out_path, ROOT)}")
    print("=" * 85)

if __name__ == '__main__':
    main()

