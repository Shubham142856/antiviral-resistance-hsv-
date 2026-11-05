"""
run_where_components_ablation.py
================================
Decomposition Ablation of the Positional Prior (`Where`):
Constraint-Only vs. Distance-Only vs. Composite Where across HSV-1 Targets.

Anchored strictly to protocol:
  docs/WHERE_COMPONENTS_ABLATION_PLAN.md (Version 1.2)

Evaluates:
  1. HSV-1 UL30-Aciclovir (Primary Arm, S_0 Pair-Weighted AUROC)
  2. HSV-1 UL30-Foscarnet (Secondary Descriptive Arm)
  3. HSV-1 UL23-Aciclovir (Comparative Non-Essential Counter-Arm)

Scores Decomposed:
  - Score 1: Evolutionary Positional Constraint Only:
      Z(Constraint) = - (mu_pLM_19 - mean) / std
  - Score 2A: 3D Active-Site Proximity Only (Linear Standardized):
      Z(Distance) = (dist - mean) / std,  Score = - Z(Distance)
  - Score 2B: 3D Active-Site Proximity Only (Rank-Normalized Percentile):
      RankNorm(Distance) = 1.0 - (rank(dist) - 1.0) / (N - 1.0)
  - Score 3A: Composite Where (Linear Standardized, Benchmark Standard):
      W_composite = (Z(Constraint) - Z(Distance)) / std(Z(Constraint) - Z(Distance))
  - Score 3B: Composite Where (Rank-Normalized Distance):
      W_rank = (Z(Constraint) + zscore(RankNorm(dist))) / std
  - Baseline: BLOSUM62 Radicalness:
      severity_blosum

Designs & Strata:
  - S0 Full Cohort Within-Study (Pair-Weighted & Study-Mean)
  - S2 Matched Quantitative EC50 Fold-Change Stratum
  - D1 Pooled Full Cohort AUROC

Resampling:
  - B = 2,000 Study-Clustered Bootstrap Resamples (seed 42)
  - Shared study resamples for exact paired deltas:
      Delta(Composite - Constraint)
      Delta(Composite - Distance)
      Delta(Distance - Constraint)

Output:
  - logs/where_components_ablation_results.json
"""

import os
import re
import json
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from scipy.stats import rankdata

warnings.filterwarnings('ignore')

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
N_BOOT = 2000
RNG_SEED = 42

def parse_mut(s):
    match = re.match(r'^([A-Z])(\d+)([A-Z])$', str(s).strip())
    if match:
        return match.group(1), int(match.group(2)), match.group(3)
    return None, None, None

def get_study_id(row):
    doi = str(row.get('ref_doi', '')).strip()
    if doi and doi not in ['nan', 'None', '']:
        return f"DOI:{doi}"
    pmid = str(row.get('PMIDs', '')).strip()
    if pmid and pmid not in ['nan', 'None', '']:
        return f"PMID:{pmid}"
    title = str(row.get('ref_title', '')).strip()
    if title and title not in ['nan', 'None', '']:
        return f"Title:{title[:40]}"
    return "Unassigned"

def compute_study_pair_stats(df_sub, score_col):
    st_groups = df_sub.groupby('study_id')
    stats = {}
    for st_id, grp in st_groups:
        res_scores = grp.loc[grp['target_binary'] == 1, score_col].values
        sens_scores = grp.loc[grp['target_binary'] == 0, score_col].values
        n_r = len(res_scores)
        n_s = len(sens_scores)
        if n_r > 0 and n_s > 0:
            n_p = n_r * n_s
            diffs = res_scores[:, None] - sens_scores[None, :]
            concordant = float(np.sum(diffs > 0) + 0.5 * np.sum(diffs == 0))
            stats[st_id] = {
                'n_res': n_r,
                'n_sens': n_s,
                'n_pairs': n_p,
                'concordant': concordant,
                'auroc': concordant / n_p
            }
    return stats

def pair_weighted_from_stats(stats_dict):
    tot_p = sum(s['n_pairs'] for s in stats_dict.values())
    tot_c = sum(s['concordant'] for s in stats_dict.values())
    return (tot_c / tot_p) if tot_p > 0 else np.nan, tot_p

def study_mean_from_stats(stats_dict):
    if not stats_dict:
        return np.nan
    return float(np.mean([s['auroc'] for s in stats_dict.values()]))

def main():
    print("=" * 85)
    print("WHERE DECOMPOSITION ABLATION: CONSTRAINT vs. DISTANCE vs. COMPOSITE")
    print("Protocol: docs/WHERE_COMPONENTS_ABLATION_PLAN.md (v1.2)")
    print("=" * 85)

    # 1. Load Primary Benchmark Cohort & Metadata
    primary_df = pd.read_csv(os.path.join(ROOT, 'data', 'processed', 'where_what_features.tsv'), sep='\t')

    # Load HerpesDRG and CMI to match Step 1 study IDs exactly
    herpes = pd.read_csv(os.path.join(ROOT, 'herpesdrg-db', 'herpesdrg-db.tsv'), sep='\t')
    hsv1_herpes = herpes[herpes['virus'] == 'HSV1'].copy()
    parsed = [parse_mut(x) for x in hsv1_herpes['aa_change']]
    hsv1_herpes['wt'] = [p[0] for p in parsed]
    hsv1_herpes['pos'] = [p[1] for p in parsed]
    hsv1_herpes['mut'] = [p[2] for p in parsed]
    hsv1_herpes['gene_std'] = hsv1_herpes['gene'].astype(str).str.strip().replace({'TK': 'UL23', 'POL': 'UL30', 'DPOL': 'UL30'})

    for col in ['strain_generation', 'test_method', 'control_strain']:
        hsv1_herpes[col] = hsv1_herpes[col].astype(str).str.strip()
        hsv1_herpes.loc[hsv1_herpes[col].isin(['nan', 'None', '']), col] = np.nan

    hdrg_cols = ['strain_generation', 'test_method', 'control_strain', 'Aciclovir', 'Foscarnet', 'ref_doi', 'ref_title']
    hsv1_herpes_dedup = hsv1_herpes.drop_duplicates(subset=['gene_std', 'pos', 'wt', 'mut']).copy()

    # Load CMI
    cmi = pd.read_csv(os.path.join(ROOT, 'data', 'processed', 'cmi_2025_all_sheets.csv'))
    cmi_hsv1 = cmi[cmi['Source_Sheet'].isin(['HSV-1 TK', 'HSV-1 Pol'])].copy()
    cmi_hsv1['gene_std'] = cmi_hsv1['Gene'].astype(str).str.strip().replace({'Thymidine Kinase': 'UL23', 'DNA Polymerase': 'UL30'})
    cmi_parsed = [parse_mut(x) for x in cmi_hsv1['AA change']]
    cmi_hsv1['wt'] = [p[0] for p in cmi_parsed]
    cmi_hsv1['pos'] = [p[1] for p in cmi_parsed]
    cmi_hsv1['mut'] = [p[2] for p in cmi_parsed]
    cmi_hsv1['drug_std'] = cmi_hsv1['Drug'].astype(str).str.strip()

    for col in ['resistance status confirmed in cell culture', 'clinical susceptibility information available', 'All genes sequenced that may confer resistance?', 'Resistance']:
        cmi_hsv1[col] = cmi_hsv1[col].astype(str).str.strip()
        cmi_hsv1.loc[cmi_hsv1[col].isin(['nan', 'None', '']), col] = np.nan

    cmi_cols = ['resistance status confirmed in cell culture', 'clinical susceptibility information available', 'All genes sequenced that may confer resistance?', 'Resistance', 'PMIDs']
    cmi_hsv1_dedup = cmi_hsv1.drop_duplicates(subset=['gene_std', 'drug_std', 'pos', 'wt', 'mut']).copy()

    # Merge
    merged = pd.merge(
        primary_df,
        hsv1_herpes_dedup[['gene_std', 'pos', 'wt', 'mut'] + hdrg_cols],
        left_on=['gene', 'pos', 'wt', 'mut'],
        right_on=['gene_std', 'pos', 'wt', 'mut'],
        how='left'
    )
    merged = pd.merge(
        merged,
        cmi_hsv1_dedup[['gene_std', 'drug_std', 'pos', 'wt', 'mut'] + cmi_cols],
        left_on=['gene', 'drug', 'pos', 'wt', 'mut'],
        right_on=['gene_std', 'drug_std', 'pos', 'wt', 'mut'],
        how='left'
    )

    def is_num_fc(row):
        val = row['Aciclovir'] if row['drug'] == 'Aciclovir' else row['Foscarnet']
        if pd.isna(val) or str(val).strip() in ['nan', 'None', '']: return False
        s = str(val).replace('>', '').replace('<', '').strip()
        try:
            float(s)
            return True
        except:
            return False

    merged['is_numeric_fc'] = [is_num_fc(r) for _, r in merged.iterrows()]
    merged['study_id'] = [get_study_id(r) for _, r in merged.iterrows()]

    # 2. Score Formulation per Gene
    for g in ['UL23', 'UL30']:
        m = (merged['gene'] == g)
        sub_g = merged.loc[m]
        
        # 1. Constraint Only: higher score = more constrained = intolerance to change
        mu = sub_g['mu_plm_19']
        z_constraint = - (mu - mu.mean()) / mu.std()
        merged.loc[m, 'score_constraint'] = z_constraint

        # 2A. Distance Only (Standardized): higher score = closer to active site
        d = sub_g['dist_to_site_angstrom']
        z_dist = (d - d.mean()) / d.std()
        merged.loc[m, 'score_dist_std'] = - z_dist

        # 2B. Distance Only (Rank-Normalized): 1.0 for closest, 0.0 for furthest
        ranks = rankdata(d, method='average')
        rank_norm = 1.0 - (ranks - 1.0) / (len(ranks) - 1.0)
        merged.loc[m, 'score_dist_ranknorm'] = rank_norm

        # 3A. Composite Where (Standardized): benchmark specification
        comp_diff = z_constraint - z_dist
        merged.loc[m, 'score_composite_std'] = comp_diff / np.std(comp_diff)

        # 3B. Composite Where (Rank-Normalized Distance)
        z_ranknorm = (rank_norm - rank_norm.mean()) / rank_norm.std()
        comp_rank = z_constraint + z_ranknorm
        merged.loc[m, 'score_composite_ranknorm'] = comp_rank / np.std(comp_rank)

        # Baseline: BLOSUM62 Radicalness
        merged.loc[m, 'score_blosum'] = sub_g['severity_blosum'].values

    # Target cohorts
    targets = [
        ('UL30-Aciclovir', 'UL30', 'Aciclovir', 'Primary Arm'),
        ('UL30-Foscarnet', 'UL30', 'Foscarnet', 'Secondary Descriptive Arm'),
        ('UL23-Aciclovir', 'UL23', 'Aciclovir', 'Comparative Counter-Arm')
    ]

    scores_to_eval = [
        ('Composite Where (Standardized)', 'score_composite_std'),
        ('Evolutionary Constraint Only (ESM-2 650M)', 'score_constraint'),
        ('Active-Site Proximity Only (Linear Std)', 'score_dist_std'),
        ('Active-Site Proximity Only (Rank-Norm)', 'score_dist_ranknorm'),
        ('Composite Where (Rank-Norm Dist)', 'score_composite_ranknorm'),
        ('BLOSUM62 Radicalness', 'score_blosum')
    ]

    results = {
        'protocol': {
            'document': 'docs/WHERE_COMPONENTS_ABLATION_PLAN.md',
            'version': '1.2',
            'status': 'COMPLETED ABLATION',
            'n_boot': N_BOOT,
            'seed': RNG_SEED
        },
        'targets': {}
    }

    rng = np.random.RandomState(RNG_SEED)

    for target_key, gene, drug, arm_role in targets:
        sub = merged[(merged['gene'] == gene) & (merged['drug'] == drug)].copy().reset_index(drop=True)
        n_total = len(sub)
        n_res = int(sub['target_binary'].sum())
        n_sens = int((sub['target_binary'] == 0).sum())

        # Pre-compute observed stats
        obs_stats = {}
        for sc_name, sc_col in scores_to_eval:
            # S0 Full Cohort
            st_dict = compute_study_pair_stats(sub, sc_col)
            s0_pw, n_pairs = pair_weighted_from_stats(st_dict)
            s0_sm = study_mean_from_stats(st_dict)

            # S2 Numeric FC Stratum
            sub_s2 = sub[sub['is_numeric_fc']].copy()
            st_dict_s2 = compute_study_pair_stats(sub_s2, sc_col)
            s2_pw, s2_pairs = pair_weighted_from_stats(st_dict_s2)
            try:
                s2_pooled = float(roc_auc_score(sub_s2['target_binary'], sub_s2[sc_col]))
            except:
                s2_pooled = np.nan

            # D1 Pooled AUROC
            try:
                d1_auc = float(roc_auc_score(sub['target_binary'], sub[sc_col]))
            except:
                d1_auc = np.nan

            obs_stats[sc_name] = {
                'col': sc_col,
                'st_dict': st_dict,
                's0_pair_weighted': s0_pw,
                's0_study_mean': s0_sm,
                'n_pairs': n_pairs,
                'm_mixed_studies': len(st_dict),
                's2_pair_weighted': s2_pw,
                's2_pooled': s2_pooled,
                's2_pairs': s2_pairs,
                's2_n_rows': len(sub_s2),
                'd1_pooled': d1_auc
            }

        # Study-clustered bootstrap with shared study resample
        st_list = list(sub['study_id'].unique())
        boot_s0 = {sc_name: [] for sc_name, _ in scores_to_eval}
        boot_deltas = {
            'Delta(Composite - Constraint)': [],
            'Delta(Composite - Distance_Std)': [],
            'Delta(Distance_Std - Constraint)': [],
            'Delta(Composite_Rank - Composite_Std)': []
        }

        for b in range(N_BOOT):
            sampled_studies = rng.choice(st_list, size=len(st_list), replace=True)
            b_vals = {}
            for sc_name, _ in scores_to_eval:
                st_dict = obs_stats[sc_name]['st_dict']
                tot_p, tot_c = 0, 0.0
                for st_id in sampled_studies:
                    if st_id in st_dict:
                        tot_p += st_dict[st_id]['n_pairs']
                        tot_c += st_dict[st_id]['concordant']
                s0_b = (tot_c / tot_p) if tot_p > 0 else np.nan
                boot_s0[sc_name].append(s0_b)
                b_vals[sc_name] = s0_b

            # Paired deltas within resample b
            w_std = b_vals['Composite Where (Standardized)']
            w_rank = b_vals['Composite Where (Rank-Norm Dist)']
            const = b_vals['Evolutionary Constraint Only (ESM-2 650M)']
            dist = b_vals['Active-Site Proximity Only (Linear Std)']

            if not np.isnan(w_std) and not np.isnan(const):
                boot_deltas['Delta(Composite - Constraint)'].append(w_std - const)
            if not np.isnan(w_std) and not np.isnan(dist):
                boot_deltas['Delta(Composite - Distance_Std)'].append(w_std - dist)
            if not np.isnan(dist) and not np.isnan(const):
                boot_deltas['Delta(Distance_Std - Constraint)'].append(dist - const)
            if not np.isnan(w_rank) and not np.isnan(w_std):
                boot_deltas['Delta(Composite_Rank - Composite_Std)'].append(w_rank - w_std)

        # Summarize target results
        target_res = {
            'gene': gene,
            'drug': drug,
            'arm_role': arm_role,
            'cohort': {
                'total_rows': n_total,
                'resistant': n_res,
                'sensitive': n_sens,
                'mixed_studies': obs_stats['Composite Where (Standardized)']['m_mixed_studies'],
                'same_study_pairs': obs_stats['Composite Where (Standardized)']['n_pairs']
            },
            'scores': {},
            'paired_deltas': {}
        }

        print(f"\nTarget: {target_key} ({arm_role})")
        print(f"  N = {n_total} ({n_res} R / {n_sens} S), Mixed Studies = {obs_stats['Composite Where (Standardized)']['m_mixed_studies']}, Pairs = {obs_stats['Composite Where (Standardized)']['n_pairs']}")
        print(f"  %-46s %8s %19s %8s %8s %8s" % ("Score", "S0_PW", "95% CI", "S0_Mean", "S2_Strat", "D1_Pooled"))
        print("  " + "-" * 105)

        for sc_name, _ in scores_to_eval:
            pw = obs_stats[sc_name]['s0_pair_weighted']
            ci = [float(np.nanpercentile(boot_s0[sc_name], 2.5)), float(np.nanpercentile(boot_s0[sc_name], 97.5))]
            sm = obs_stats[sc_name]['s0_study_mean']
            s2 = obs_stats[sc_name]['s2_pair_weighted']
            s2_pooled = obs_stats[sc_name]['s2_pooled']
            d1 = obs_stats[sc_name]['d1_pooled']

            ci_str = f"[{ci[0]:.4f}, {ci[1]:.4f}]"
            s2_str = f"{s2_pooled:.4f}" if not np.isnan(s2_pooled) else "n/a"
            print(f"  %-46s %8.4f %19s %8.4f %8s %8.4f" % (sc_name, pw, ci_str, sm, s2_str, d1))

            target_res['scores'][sc_name] = {
                's0_pair_weighted': round(pw, 4),
                's0_pair_weighted_ci_95': [round(ci[0], 4), round(ci[1], 4)],
                's0_study_mean': round(sm, 4),
                's2_stratum_pair_weighted': round(s2, 4) if not np.isnan(s2) else None,
                's2_stratum_pooled_auc': round(s2_pooled, 4) if not np.isnan(s2_pooled) else None,
                'd1_pooled': round(d1, 4)
            }

        print("  Paired Bootstrap Deltas (Within Shared Study Resamples):")
        for d_name, d_arr in boot_deltas.items():
            obs_delta = float(np.nanmean(d_arr))
            d_ci = [float(np.nanpercentile(d_arr, 2.5)), float(np.nanpercentile(d_arr, 97.5))]
            p_val = float(np.mean(np.array(d_arr) <= 0)) if obs_delta > 0 else float(np.mean(np.array(d_arr) >= 0))
            print(f"    {d_name:<38} : {obs_delta:+.4f} [{d_ci[0]:+.4f}, {d_ci[1]:+.4f}] (p = {p_val:.4f})")

            target_res['paired_deltas'][d_name] = {
                'observed_delta': round(obs_delta, 4),
                'ci_95': [round(d_ci[0], 4), round(d_ci[1], 4)],
                'p_value': round(p_val, 4)
            }

        results['targets'][target_key] = target_res

    # Save to logs
    out_path = os.path.join(ROOT, 'logs', 'where_components_ablation_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)

    print("\n" + "=" * 85)
    print(f"Wrote complete ablation results to: {os.path.relpath(out_path, ROOT)}")
    print("=" * 85)

if __name__ == '__main__':
    main()
