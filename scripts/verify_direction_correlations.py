import pathlib
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

def main():
    ds_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/primary_position_folds.tsv')
    cache_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/esm2_pos_cache.json')
    
    df = pd.read_csv(ds_path, sep='\t')
    with open(cache_path, 'r', encoding='utf-8') as f:
        raw_cache = json.load(f)
    pos_llr_cache = { (k.split(':')[0], int(k.split(':')[1])): v for k, v in raw_cache.items() }
    
    llrs = []
    for idx, r in df.iterrows():
        g, p, mut = r['gene'], int(r['pos']), r['mut']
        scores = pos_llr_cache.get((g, p), {})
        llrs.append(scores.get(mut, 0.0))
    df['esm2_llr'] = llrs
    
    pairs = [
        ('UL23', 'Aciclovir', 'UL23-Aciclovir'),
        ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
        ('UL30', 'Foscarnet', 'UL30-Foscarnet')
    ]
    
    results = {}
    lines = [
        "=== ESM-2 (35M) DELTA-LLR BIOLOGICAL DIRECTION CORRELATION AUDIT ===",
        "Delta-LLR Definition: log P(mut | x_{\\i}) - log P(wt | x_{\\i})",
        "Hypothesis 1: Negative correlation (Delta-LLR < 0, deleterious -> resistant)",
        "Hypothesis 2: Positive correlation (Delta-LLR > 0, stabilizing/fit -> resistant)",
        ""
    ]
    
    for g, d, p_key in pairs:
        sub = df[(df['gene'] == g) & (df['drug'] == d)].copy().reset_index(drop=True)
        y = sub['target_binary'].values.astype(int)
        scores = sub['esm2_llr'].values
        
        # Split into sensitive and resistant
        sens_llrs = scores[y == 0]
        res_llrs = scores[y == 1]
        
        # Correlations
        spearman_rho, sp_p = stats.spearmanr(scores, y)
        kendall_tau, kd_p = stats.kendalltau(scores, y)
        pt_biserial, pb_p = stats.pointbiserialr(y, scores)
        
        # Mann-Whitney U test (testing if resistant has lower scores than sensitive)
        mwu_stat, mwu_p = stats.mannwhitneyu(res_llrs, sens_llrs, alternative='less')
        
        auc_raw = roc_auc_score(y, scores)
        auc_inv = roc_auc_score(y, -scores)
        
        is_negative = bool(spearman_rho < 0 and mwu_p < 0.05)
        is_positive = bool(spearman_rho > 0 and mwu_p > 0.95)
        
        if is_negative:
            verdict = "Consistently Deleterious -> Resistant (Negative Correlation)"
        elif is_positive:
            verdict = "Consistently Fit/Stabilizing -> Resistant (Positive Correlation)"
        else:
            verdict = "Mixed / Weak Correlation"
            
        lines.append(f"=== {p_key} (N={len(sub)}: {len(res_llrs)} Resistant, {len(sens_llrs)} Sensitive) ===")
        lines.append(f"  Mean Delta-LLR Sensitive: {np.mean(sens_llrs):.4f} +/- {np.std(sens_llrs):.4f}")
        lines.append(f"  Mean Delta-LLR Resistant: {np.mean(res_llrs):.4f} +/- {np.std(res_llrs):.4f}")
        lines.append(f"  Spearman rho:             {spearman_rho:+.4f} (p = {sp_p:.4e})")
        lines.append(f"  Kendall tau:              {kendall_tau:+.4f} (p = {kd_p:.4e})")
        lines.append(f"  Point-Biserial r:         {pt_biserial:+.4f} (p = {pb_p:.4e})")
        lines.append(f"  Mann-Whitney U (Res<Sens): p = {mwu_p:.4e}")
        lines.append(f"  AUROC (Raw Delta-LLR):    {auc_raw:.4f}")
        lines.append(f"  AUROC (-Delta-LLR):       {auc_inv:.4f}")
        lines.append(f"  Direction Verdict:        {verdict}")
        lines.append("")
        
        results[p_key] = {
            'mean_llr_sensitive': round(float(np.mean(sens_llrs)), 4),
            'std_llr_sensitive': round(float(np.std(sens_llrs)), 4),
            'mean_llr_resistant': round(float(np.mean(res_llrs)), 4),
            'std_llr_resistant': round(float(np.std(res_llrs)), 4),
            'spearman_rho': round(float(spearman_rho), 4),
            'spearman_p': float(sp_p),
            'kendall_tau': round(float(kendall_tau), 4),
            'kendall_p': float(kd_p),
            'point_biserial_r': round(float(pt_biserial), 4),
            'point_biserial_p': float(pb_p),
            'mwu_res_less_sens_p': float(mwu_p),
            'auroc_raw': round(float(auc_raw), 4),
            'auroc_inverted': round(float(auc_inv), 4),
            'verdict': verdict
        }
        
    out_json = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/direction_correlations_results.json')
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
        
    out_log = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/direction_correlations.log')
    with open(out_log, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines))
        
    print("\n".join(lines))

if __name__ == '__main__':
    main()

