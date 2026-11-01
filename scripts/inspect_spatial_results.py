import pathlib
import json

with open(str(pathlib.Path(__file__).resolve().parent.parent / 'logs/spatial_3d_holdout_results.json')) as f:
    d = json.load(f)

for pair in d['results_by_pair']:
    print(f"\n========================================================")
    print(f"=== {pair} ===")
    print(f"========================================================")
    for cohort in ['full_cohort', 'resolved_core']:
        print(f"\n--- Cohort: {cohort} ---")
        c_data = d['results_by_pair'][pair][cohort]
        for m in c_data:
            b1d = c_data[m]['1d_block']
            s3d = c_data[m]['3d_spatial']
            delta = c_data[m]['delta_3d_minus_1d']
            print(f"{m:<38} | 1D: {b1d['mean_auc']:.4f} (Pooled {b1d['pooled_auc']:.4f} [{b1d['pooled_95ci'][0]:.4f}, {b1d['pooled_95ci'][1]:.4f}]) | 3D: {s3d['mean_auc']:.4f} (Pooled {s3d['pooled_auc']:.4f} [{s3d['pooled_95ci'][0]:.4f}, {s3d['pooled_95ci'][1]:.4f}]) | Delta: {delta['mean_diff']:+.4f}")

