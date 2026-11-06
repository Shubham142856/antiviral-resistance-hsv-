"""
scripts/compute_pssm_esm2_paired_contrast.py
============================================
Computes the paired within-study (D3) contrast between multi-herpesvirus PSSM
and ESM-2 650M under a shared study-clustered bootstrap (B=2,000, seed=42).

Addresses the difference-of-significance fallacy (Gelman & Stern, 2006):
Observing that PSSM's 95% CI strictly excludes 0.50 while ESM-2's CI spans 0.50
(on UL30-FOS) does not establish a statistically significant difference between them.
This script directly tests the paired contrast Delta = S0(PSSM) - S0(ESM-2)
on identical resampled study clusters.
"""

import os
import re
import json
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT_PATH = os.path.join(ROOT, "logs", "pssm_vs_esm2_d3_paired_contrast.json")

def parse_mut(s):
    match = re.match(r'^([A-Z])(\d+)([A-Z])$', str(s).strip())
    if match:
        return match.group(1), int(match.group(2)), match.group(3)
    return None, None, None

def get_study_id(row):
    doi = str(row['ref_doi']).strip()
    if doi and doi not in ['nan', 'None', '']:
        return f"DOI:{doi}"
    pmid = str(row['PMIDs']).strip()
    if pmid and pmid not in ['nan', 'None', '']:
        return f"PMID:{pmid}"
    title = str(row['ref_title']).strip()
    if title and title not in ['nan', 'None', '']:
        return f"Title:{title[:40]}"
    return "Unassigned"

def load_data():
    df = pd.read_csv(os.path.join(ROOT, 'data/processed/where_what_features.tsv'), sep='\t')
    spatial_df = pd.read_csv(os.path.join(ROOT, 'data/processed/primary_3d_spatial_folds.tsv'), sep='\t')
    cmi_all = pd.read_csv(os.path.join(ROOT, 'data/processed/cmi_2025_all_sheets.csv'))
    herpes = pd.read_csv(os.path.join(ROOT, 'herpesdrg-db/herpesdrg-db.tsv'), sep='\t')
    
    cols_to_merge = ['score_esm2_650m', 'score_pssm']
    df = pd.merge(
        df,
        spatial_df[['gene', 'pos', 'wt', 'mut', 'drug'] + cols_to_merge],
        on=['gene', 'pos', 'wt', 'mut', 'drug'],
        how='left'
    )
    
    hsv1_herpes = herpes[herpes['virus'] == 'HSV1'].copy()
    parsed = [parse_mut(x) for x in hsv1_herpes['aa_change']]
    hsv1_herpes['wt'] = [p[0] for p in parsed]
    hsv1_herpes['pos'] = [p[1] for p in parsed]
    hsv1_herpes['mut'] = [p[2] for p in parsed]
    hsv1_herpes['gene_std'] = hsv1_herpes['gene'].astype(str).str.strip().replace({'TK': 'UL23', 'POL': 'UL30', 'DPOL': 'UL30'})
    hsv1_herpes_dedup = hsv1_herpes.drop_duplicates(subset=['gene_std', 'pos', 'wt', 'mut']).copy()
    
    df = pd.merge(
        df,
        hsv1_herpes_dedup[['gene_std', 'pos', 'wt', 'mut', 'Aciclovir', 'Foscarnet', 'ref_doi', 'ref_title']],
        left_on=['gene', 'pos', 'wt', 'mut'],
        right_on=['gene_std', 'pos', 'wt', 'mut'],
        how='left'
    )
    
    cmi_hsv1 = cmi_all[cmi_all['Source_Sheet'].isin(['HSV-1 TK', 'HSV-1 Pol'])].copy()
    cmi_hsv1['gene_std'] = cmi_hsv1['Gene'].astype(str).str.strip().replace({'Thymidine Kinase': 'UL23', 'DNA Polymerase': 'UL30'})
    c_parsed = [parse_mut(x) for x in cmi_hsv1['AA change']]
    cmi_hsv1['wt'] = [p[0] for p in c_parsed]
    cmi_hsv1['pos'] = [p[1] for p in c_parsed]
    cmi_hsv1['mut'] = [p[2] for p in c_parsed]
    cmi_hsv1['drug_std'] = cmi_hsv1['Drug'].astype(str).str.strip()
    cmi_hsv1_dedup = cmi_hsv1.drop_duplicates(subset=['gene_std', 'drug_std', 'pos', 'wt', 'mut']).copy()
    
    df = pd.merge(
        df,
        cmi_hsv1_dedup[['gene_std', 'drug_std', 'pos', 'wt', 'mut', 'PMIDs']],
        left_on=['gene', 'drug', 'pos', 'wt', 'mut'],
        right_on=['gene_std', 'drug_std', 'pos', 'wt', 'mut'],
        how='left'
    )
    
    df['study_id'] = [get_study_id(r) for _, r in df.iterrows()]
    return df

def main():
    df = load_data()
    
    TARGETS = [
        ('UL30', 'Aciclovir', 'UL30-Aciclovir'),
        ('UL30', 'Foscarnet', 'UL30-Foscarnet'),
        ('UL23', 'Aciclovir', 'UL23-Aciclovir')
    ]
    
    SCORE_COLS = {
        'PSSM': 'score_pssm',
        'ESM-2 650M': 'score_esm2_650m'
    }
    
    # 1. Identify all mixed studies across targets
    target_s0_data = {}
    all_mixed_studies = set()
    
    for g, d, t_key in TARGETS:
        sub = df[(df['gene'] == g) & (df['drug'] == d)].copy()
        st_groups = sub.groupby('study_id')
        
        t_scores_data = {}
        for s_label, s_col in SCORE_COLS.items():
            s_stats = {}
            for st_id, grp in st_groups:
                rv = grp.loc[grp['target_binary'] == 1, s_col].values.astype(float)
                sv = grp.loc[grp['target_binary'] == 0, s_col].values.astype(float)
                nr = len(rv)
                ns = len(sv)
                if nr > 0 and ns > 0:
                    n_p = nr * ns
                    diffs = rv[:, None] - sv[None, :]
                    conc = float(np.sum(diffs > 0) + 0.5 * np.sum(diffs == 0))
                    s_stats[st_id] = {'pairs': n_p, 'conc': conc, 'auc': conc / n_p}
                    all_mixed_studies.add(st_id)
            t_scores_data[s_label] = s_stats
        target_s0_data[t_key] = t_scores_data
        
    study_universe = sorted(list(all_mixed_studies))
    n_studies = len(study_universe)
    study_to_idx = {st: i for i, st in enumerate(study_universe)}
    print(f"Total mixed studies in unified universe: {n_studies}")
    
    # Precompute pair and conc arrays
    pair_arrays = {}
    conc_arrays = {s_label: {} for s_label in SCORE_COLS}
    
    for g, d, t_key in TARGETS:
        p_arr = np.zeros(n_studies)
        for st_id, s_data in target_s0_data[t_key]['PSSM'].items():
            p_arr[study_to_idx[st_id]] = s_data['pairs']
        pair_arrays[t_key] = p_arr
        
        for s_label in SCORE_COLS:
            c_arr = np.zeros(n_studies)
            for st_id, s_data in target_s0_data[t_key][s_label].items():
                c_arr[study_to_idx[st_id]] = s_data['conc']
            conc_arrays[s_label][t_key] = c_arr
            
    # Point estimates
    point_estimates = {}
    for g, d, t_key in TARGETS:
        tot_pairs = np.sum(pair_arrays[t_key])
        pssm_conc = np.sum(conc_arrays['PSSM'][t_key])
        esm_conc = np.sum(conc_arrays['ESM-2 650M'][t_key])
        
        pssm_s0 = float(pssm_conc / tot_pairs) if tot_pairs > 0 else np.nan
        esm_s0 = float(esm_conc / tot_pairs) if tot_pairs > 0 else np.nan
        delta_obs = pssm_s0 - esm_s0
        
        point_estimates[t_key] = {
            'total_pairs': int(tot_pairs),
            'pssm_s0': pssm_s0,
            'esm2_650m_s0': esm_s0,
            'observed_delta_pssm_minus_esm2': delta_obs
        }
        
    # Shared study-clustered bootstrap (B=2000, seed=42)
    np.random.seed(42)
    B = 2000
    boot_deltas = {t_key: [] for _, _, t_key in TARGETS}
    boot_pssm = {t_key: [] for _, _, t_key in TARGETS}
    boot_esm = {t_key: [] for _, _, t_key in TARGETS}
    
    for b in range(B):
        resampled_idx = np.random.choice(n_studies, size=n_studies, replace=True)
        counts = np.bincount(resampled_idx, minlength=n_studies)
        
        for g, d, t_key in TARGETS:
            p_tot = np.dot(counts, pair_arrays[t_key])
            if p_tot > 0:
                p_s0 = float(np.dot(counts, conc_arrays['PSSM'][t_key]) / p_tot)
                e_s0 = float(np.dot(counts, conc_arrays['ESM-2 650M'][t_key]) / p_tot)
                d_b = p_s0 - e_s0
                boot_deltas[t_key].append(d_b)
                boot_pssm[t_key].append(p_s0)
                boot_esm[t_key].append(e_s0)
            else:
                boot_deltas[t_key].append(np.nan)
                boot_pssm[t_key].append(np.nan)
                boot_esm[t_key].append(np.nan)
                
    results = {
        'meta': {
            'description': 'Paired D3 within-study contrast: PSSM minus ESM-2 650M',
            'n_boot': B,
            'seed': 42,
            'n_mixed_studies': n_studies,
            'status': 'EXPLORATORY'
        },
        'targets': {}
    }
    
    print("\n" + "="*80)
    print("PAIRED D3 CONTRAST RESULTS: PSSM minus ESM-2 650M (B=2,000, seed=42)")
    print("="*80)
    
    for g, d, t_key in TARGETS:
        clean_d = [x for x in boot_deltas[t_key] if not np.isnan(x)]
        ci_95 = [float(np.percentile(clean_d, 2.5)), float(np.percentile(clean_d, 97.5))]
        
        obs_d = point_estimates[t_key]['observed_delta_pssm_minus_esm2']
        # Two-sided empirical p-value against null of delta = 0
        p_two_sided = float(2.0 * min(np.mean(np.array(clean_d) <= 0), np.mean(np.array(clean_d) >= 0)))
        p_two_sided = min(p_two_sided, 1.0)
        p_str = f"p < {1.0/B:.4f}" if p_two_sided < 1.0/B else f"p = {p_two_sided:.4f}"
        
        excludes_zero = bool(ci_95[0] > 0 or ci_95[1] < 0)
        
        target_res = {
            'target': t_key,
            'total_pairs': point_estimates[t_key]['total_pairs'],
            'pssm_d3_s0': point_estimates[t_key]['pssm_s0'],
            'esm2_650m_d3_s0': point_estimates[t_key]['esm2_650m_s0'],
            'observed_delta': round(obs_d, 4),
            'ci_95': [round(ci_95[0], 4), round(ci_95[1], 4)],
            'p_value_two_sided': p_two_sided,
            'p_value_formatted': p_str,
            'excludes_zero': excludes_zero,
            'verdict': 'STATISTICALLY SIGNIFICANT' if excludes_zero else 'INCONCLUSIVE (CI spans zero)'
        }
        results['targets'][t_key] = target_res
        
        print(f"\nTarget: {t_key}")
        print(f"  Total within-study pairs: {target_res['total_pairs']}")
        print(f"  PSSM D3 S0:       {target_res['pssm_d3_s0']:.4f}")
        print(f"  ESM-2 650M D3 S0: {target_res['esm2_650m_d3_s0']:.4f}")
        print(f"  Paired Delta (PSSM - ESM-2): {target_res['observed_delta']:+.4f} [{target_res['ci_95'][0]:+.4f}, {target_res['ci_95'][1]:+.4f}]")
        print(f"  Two-sided p-value: {p_str}")
        print(f"  Verdict: {target_res['verdict']}")
        
    with open(OUT_PATH, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved results to {OUT_PATH}")

if __name__ == '__main__':
    main()
