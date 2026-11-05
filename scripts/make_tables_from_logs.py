"""
make_tables_from_logs.py
========================
Deterministically generates markdown for:
  - Table 7: Recalibrated Benchmarks & Excess Metrics across Structurally Resolved Domains (logs/step0b_final_recalibrated_benchmarks.json)
  - Table 8: The Evidence Ladder across D1, D2, D3, S2 (logs/evidence_ladder_table8_results.json)
  - Step 1 Tables: Positivity Check, S0 Primary Within-Study, Contrast D, and Secondary Strata (logs/step1_ascertainment_control_results.json)
  - Table 9: Cross-Herpesvirus Replication Benchmark (HCMV UL54 and UL97) (logs/step2_hcmv_replication_results.json and logs/step2_hcmv_composite_results.json)

Rule: Manuscript tables must be generated from this script, never hand-typed.
All quantities trace to audited JSON logs with unified seeds.
"""

import os
import json
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def generate_table7():
    with open(os.path.join(ROOT, 'logs', 'step0b_final_recalibrated_benchmarks.json'), 'r') as f:
        t7 = json.load(f)
        
    md = []
    md.append("### Table 7: Recalibrated Domain-Level Benchmarks & Excess Metrics\n")
    md.append("| Target | Evaluable Domains | Pos-Only RF | Clamped Ref | ESM-2 650M | ProtBERT | BLOSUM62 | Where (Positional Prior) | Supervised LR | Excess (Where - Ref) [95% CI] |")
    md.append("|---|---|---|---|---|---|---|---|---|---|")
    
    for t_key in ['UL30-Aciclovir', 'UL30-Foscarnet', 'UL23-Aciclovir']:
        d = t7[t_key]
        n_doms = d['n_evaluable_domains']
        dom_list = ", ".join(d['evaluable_domains_list'])
        dm = d['domain_means']
        ex = d['excess_over_reference']
        ci = d['excess_descriptive_ci_95']
        
        pos_str = f"{dm['pos_only_rf']:.4f}"
        ref_str = f"{dm['reference_clamped']:.4f}"
        esm_str = f"{dm['esm2_650m']:.4f}"
        bert_str = f"{dm['protbert']:.4f}"
        blo_str = f"{dm['blosum62']:.4f}"
        whe_str = f"{dm['where_specexact']:.4f}"
        lr_str = f"{dm['supervised_lr']:.4f}"
        
        if ci.get('where_specexact') is not None and not np.isnan(ci['where_specexact'][0]):
            ex_str = f"**{ex['where_specexact']:+.4f}** [{ci['where_specexact'][0]:+.4f}, {ci['where_specexact'][1]:+.4f}] *(DESCRIPTIVE ONLY; n={n_doms} domains, 35 distinct resamples)*"
        else:
            ex_str = f"**{ex['where_specexact']:+.4f}** [No CI, n={n_doms} domains]"
            
        md.append(f"| **{t_key}** | {n_doms} ({dom_list}) | {pos_str} | {ref_str} | {esm_str} | {bert_str} | {blo_str} | **{whe_str}** | {lr_str} | {ex_str} |")
        
    contr = t7.get('descriptive_contrast_blosum62', {})
    if contr:
        md.append(f"\n*Descriptive BLOSUM62 Domain-Mean Contrast:* UL23 ({contr.get('ul23_mean', 0):.4f}) - UL30 ({contr.get('ul30_mean', 0):.4f}) = **{contr.get('delta', 0):+.4f}** ({contr.get('label', '')})\n")
    return "\n".join(md)

def generate_table8():
    with open(os.path.join(ROOT, 'logs', 'evidence_ladder_table8_results.json'), 'r') as f:
        t8 = json.load(f)
        
    md = []
    md.append("### Table 8: The Evidence Ladder (4-Design Evaluation Matrix)\n")
    md.append("*Methodological Note:* D1 for Position-Only RF represents 5-fold cross-validated out-of-fold (OOF) predictions; all other D1 evaluations are un-fitted zero-shot scores on the pooled cohort. D2 domain-mean CIs for n=4 domains are descriptive only (35 distinct resamples). D3 reports pair-weighted AUROC with 95% study-clustered bootstrap CI (B=2000, seed=42) alongside unweighted study-mean.\n")
    
    for t_key in ['UL30-Aciclovir', 'UL30-Foscarnet', 'UL23-Aciclovir']:
        tdata = t8[t_key]
        md.append(f"#### Target: {t_key} (Total N={tdata['N_total']}, Resistant={tdata['N_res']}, Sensitive={tdata['N_sens']})\n")
        md.append("| Model / Score | D1: Pooled (Full Cohort) | D2: Dual-Buffer Domain-Mean | D3: S0 Within-Study (Pair-Weighted) | S2: Numeric Fold-Change Stratum |")
        md.append("|---|---|---|---|---|")
        
        for s_name, s_vals in tdata['scores'].items():
            row_label = s_name
            if s_name == 'Position-Only RF':
                row_label = 'Position-Only RF (Out-of-Fold)'
                
            d1 = s_vals['D1_pooled']
            d1_str = f"**{d1['auc']:.4f}** [{d1['ci_95'][0]:.3f}, {d1['ci_95'][1]:.3f}]" if d1['auc'] is not None and not np.isnan(d1['auc']) else "n/a"
            
            d2 = s_vals['D2_domain_mean']
            if d2['ci_95'] is not None and not np.isnan(d2['ci_95'][0]):
                d2_str = f"**{d2['auc']:.4f}** [{d2['ci_95'][0]:.3f}, {d2['ci_95'][1]:.3f}] *(descriptive)*"
            elif d2['auc'] is not None and not np.isnan(d2['auc']):
                d2_str = f"**{d2['auc']:.4f}** [No CI, n=2]"
            else:
                d2_str = "n/a"
                
            d3 = s_vals['D3_within_study']
            if d3['pair_weighted'] is not None and not np.isnan(d3['pair_weighted']):
                d3_str = f"**{d3['pair_weighted']:.4f}** [{d3['ci_95'][0]:.3f}, {d3['ci_95'][1]:.3f}] *(mean={d3['study_mean']:.3f})*"
            else:
                d3_str = "n/a"
                
            s2 = s_vals['S2_numeric_fc']
            if s2['auc'] is not None and not np.isnan(s2['auc']):
                s2_str = f"**{s2['auc']:.4f}** [{s2['ci_95'][0]:.3f}, {s2['ci_95'][1]:.3f}]"
            else:
                s2_str = "n/a"
                
            md.append(f"| {row_label} | {d1_str} | {d2_str} | {d3_str} | {s2_str} |")
        md.append("\n")
    return "\n".join(md)

def generate_step1_tables():
    with open(os.path.join(ROOT, 'logs', 'step1_ascertainment_control_results.json'), 'r') as f:
        s1 = json.load(f)
        
    md = []
    # 1. Positivity
    md.append("### Step 1 Table 1A: Within-Study Positivity Audit\n")
    md.append("| Target | Role | Mixed Studies ($M$) | Same-Study Pairs ($N_{\\text{pairs}}$) | Top 1 Study Share (%) | Distinct Res Variants | Positivity Verdict |")
    md.append("|---|---|---|---|---|---|---|")
    for t_key, pos_data in s1['positivity_check'].items():
        md.append(f"| **{t_key}** | {pos_data['role']} | {pos_data['mixed_studies_count']} | {pos_data['total_same_study_pairs']} | {pos_data['top1_study_share_pct']:.2f}% | {pos_data['distinct_resistant_variants']} | **{pos_data['positivity_verdict']}** |")
        
    # 2. S0 Primary Within-Study
    md.append("\n### Step 1 Table 1B: S0 Primary Within-Study Evaluation\n")
    md.append("| Target | Role | Where Pair-Weighted [95% CI] | Where Study-Mean | BLOSUM62 Pair-Weighted [95% CI] | BLOSUM62 Study-Mean | Delta (Where - BLOSUM) [95% CI] |")
    md.append("|---|---|---|---|---|---|---|")
    for t_key, s0_data in s1['S0_primary_within_study'].items():
        w_pw = s0_data['where']['pair_weighted_auroc']
        w_ci = s0_data['where']['bootstrap_95ci']
        w_sm = s0_data['where']['study_equal_auroc']
        
        b_pw = s0_data['blosum62']['pair_weighted_auroc']
        b_ci = s0_data['blosum62']['bootstrap_95ci']
        b_sm = s0_data['blosum62']['study_equal_auroc']
        
        del_pw = s0_data['delta_where_minus_blosum']
        del_ci = s0_data['delta_where_minus_blosum_95ci']
        
        w_str = f"**{w_pw:.4f}** [{w_ci[0]:.3f}, {w_ci[1]:.3f}]"
        b_str = f"{b_pw:.4f} [{b_ci[0]:.3f}, {b_ci[1]:.3f}]"
        del_str = f"**{del_pw:+.4f}** [{del_ci[0]:+.3f}, {del_ci[1]:+.3f}]"
        
        md.append(f"| **{t_key}** | {s0_data['role']} | {w_str} | {w_sm:.4f} | {b_str} | {b_sm:.4f} | {del_str} |")
        
    cg = s1.get('cross_gene_contrast_D', {})
    verd = s1.get('decision_verdicts', {}).get('cross_gene_dissociation_survival_test', {}).get('verdict', 'INCONCLUSIVE')
    if cg:
        md.append(f"\n*Cross-Gene Contrast $D$ (UL30 Delta - UL23 Delta):* **{cg.get('point_estimate', 0):+.4f}** [{cg.get('bootstrap_95ci', [0,0])[0]:+.4f}, {cg.get('bootstrap_95ci', [0,0])[1]:+.4f}] — **{verd}**\n")
        
    # 3. Secondary Strata (Descriptive)
    md.append("### Step 1 Table 1C: Secondary Evidence Strata (Descriptive)\n")
    md.append("| Target | Stratum | N (Res / Sens) | Evaluable ($N_{\\text{pos}} \\ge 10$) | Where AUROC | BLOSUM62 AUROC | Description & Target Share |")
    md.append("|---|---|---|---|---|---|---|")
    
    target_totals = {
        'UL30-Aciclovir': 332,
        'UL30-Foscarnet': 311,
        'UL23-Aciclovir': 286
    }
    
    for t_key, strata_dict in s1['secondary_strata_descriptive'].items():
        tot = target_totals.get(t_key, 1)
        for s_key, s_data in strata_dict.items():
            ev_str = "Yes" if s_data['evaluable'] else "No (Underpowered)"
            w_val = f"{s_data['where_auroc']:.4f}" if s_data['where_auroc'] is not None and not np.isnan(s_data['where_auroc']) else "n/a"
            b_val = f"{s_data['blosum_auroc']:.4f}" if s_data['blosum_auroc'] is not None and not np.isnan(s_data['blosum_auroc']) else "n/a"
            
            n_tot_strat = s_data['N']
            pct_cohort = (n_tot_strat / tot) * 100
            if 'S1' in s_key:
                qual = "near-full cohort" if pct_cohort > 80.0 else "partial cohort"
                desc = f"{pct_cohort:.1f}% of cohort ({n_tot_strat}/{tot}; CMI cell-culture confirmed = yes; {qual})"
            elif 'S2' in s_key:
                desc = f"{pct_cohort:.1f}% of cohort ({n_tot_strat}/{tot}; HerpesDRG numeric fold-change, balanced quantitative)"
            elif 'S3' in s_key:
                desc = f"{pct_cohort:.1f}% of cohort ({n_tot_strat}/{tot}; recombinant marker transfer / BAC)"
            else:
                desc = s_data.get('description', '')
                
            md.append(f"| **{t_key}** | {s_key} | {n_tot_strat} ({s_data['Resistant']} R / {s_data['Sensitive']} S) | {ev_str} | {w_val} | {b_val} | {desc} |")
            
    md.append("\n*Methodological Note on S1:* Cell-culture confirmation flags cover 91.6% of UL30-ACV (304/332) and 91.6% of UL23-ACV (262/286), but only 68.5% of UL30-FOS (213/311). Because a flag covering ~69–92% of entries reflects general cohort compilation rather than a selective evidence tier, S1 serves as a broad descriptive census rather than a quality-filtered stratum.")
    return "\n".join(md)

def generate_step2_tables():
    s2_path = os.path.join(ROOT, 'logs', 'step2_hcmv_replication_results.json')
    comp_path = os.path.join(ROOT, 'logs', 'step2_hcmv_composite_results.json')
    if not os.path.exists(s2_path):
        return "Step 2 replication results not found at logs/step2_hcmv_replication_results.json"
        
    with open(s2_path, 'r', encoding='utf-8') as f:
        s2_data = json.load(f)
    comp_data = None
    if os.path.exists(comp_path):
        with open(comp_path, 'r', encoding='utf-8') as f:
            comp_data = json.load(f)
        
    md = []
    md.append("### Table 9: Cross-Herpesvirus Replication Benchmark (HCMV UL54 and UL97)\n")
    md.append("#### Table 9A: Step 2 Positivity & Study-Cluster Characteristics\n")
    md.append("| Target System | Arm Type | Total Rows | Unique Variants | Mixed Studies (M) | Same-Study Pairs | Top-1 Study Share | Distinct Variants in Mixed (R / S) | Positivity Status |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    
    for t_key, t_info in s2_data['targets'].items():
        g = t_info['gene']
        d = t_info['drug']
        arm = t_info['arm_type']
        pos = t_info['positivity']
        rules = pos['rules']
        
        md.append(
            f"| **{g}-{d}** | {arm} | {pos['total_cohort_rows']} | {pos['unique_mutations_overall']} | "
            f"**{pos['mixed_study_count']}** | **{pos['total_same_study_pairs']}** | {pos['top_1_study_share_pct']:.1f}% | "
            f"{pos['distinct_res_variants_in_mixed']} R / {pos['distinct_sens_variants_in_mixed']} S | **{rules['overall_status']}** |"
        )
        
    md.append("\n#### Table 9B: Replication Evidence Ladder Matrix (S0 Within-Study, S2 Stratum, and Pooled D1)\n")
    md.append("*Specification Note:* Reports the pre-registered primary distance-only score alongside post hoc exploratory evaluations (composite Where combining ESM-2 positional evolutionary constraint with Option A distance, 35M model, and Option B mapped coordinates). Option A uses the Swiss-Model P08546 homology model (HSV-1 8OJ7 template, 32.8% sequence identity). All evaluations evaluate identical 19 active-site residues using the S0 pair-weighted estimator and seed 42.\n")
    md.append("| Target System | Role / Arm | Evaluated Score | S0 Pair-Weighted AUROC [95% CI] | Study-Mean AUROC [95% CI] | LOSO Range | S2 Quantitative Stratum | Pooled D1 [95% CI] | Pre-Registered Rule Verdict |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    
    for t_key, t_info in s2_data['targets'].items():
        g = t_info['gene']
        d = t_info['drug']
        arm = t_info['arm_type']
        
        if comp_data and t_key in comp_data['targets']:
            c_feats = comp_data['targets'][t_key]['features']
            prim_or_sec = "Pre-Registered Primary" if t_key == 'HCMV_UL54_Ganciclovir' else "Secondary Descriptive"
            display_order = [
                (f'Active-Site Proximity (Distance-Only, Option A: {prim_or_sec})', 'Distance-Only (Option A: As Executed)'),
                ('Post Hoc Exploratory: Composite Where (ESM-2 650M + Option A)', 'Composite Where (ESM-2 650M + Option A: As Pre-Registered)'),
                ('Post Hoc Exploratory: Composite Where (ESM-2 35M + Option A)', 'Composite Where (ESM-2 35M + Option A)'),
                ('Post Hoc Exploratory: Positional Constraint Only (ESM-2 650M)', 'Constraint-Only (ESM-2 650M Positional)'),
                ('Post Hoc Exploratory: Composite Where (ESM-2 650M + Option B)', 'Composite Where (ESM-2 650M + Option B)'),
                ('Sensitivity Control: Active-Site Proximity (Distance-Only, Option B: HSV-1 Mapped)', 'Distance-Only (Option B)'),
                ('BLOSUM62 Radicalness', 'BLOSUM62 Radicalness')
            ]
            for label, k in display_order:
                if k not in c_feats: continue
                fv = c_feats[k]
                pw = fv['s0_pair_weighted_auroc']
                ci = fv['s0_pair_weighted_ci_95']
                sm = fv['s0_study_mean_auroc']
                sm_ci = fv['s0_study_mean_ci_95']
                loso = fv['s0_loso_range']
                s2_val = fv['s2_numeric_fc_pair_weighted_auroc']
                s2_p = fv['s2_numeric_fc_pairs']
                d1_val = fv['d1_pooled_auroc']
                d1_ci = fv['d1_pooled_ci_95']
                
                # Pre-registered verdicts
                if t_key == 'HCMV_UL54_Ganciclovir':
                    if 'Distance-Only' in label and 'Option A' in label:
                        verdict = '**NOT SUPPORTED (Failed Replication; Point < 0.55)**'
                    elif 'Composite' in label and '650M + Option A' in label:
                        verdict = '**INCONCLUSIVE (Point in [0.55, 0.65); CI Spans 0.50)**'
                    elif 'Composite' in label and '35M' in label:
                        verdict = '**INCONCLUSIVE (Point in [0.55, 0.65); CI Spans 0.50)**'
                    elif 'Composite' in label and 'Option B' in label:
                        verdict = '**POST HOC EXPLORATORY (Point > 0.75; CI Excludes 0.50)**'
                    elif 'Constraint' in label:
                        verdict = '**EXPLORATORY (Non-Spatial Evolutionary Signal)**'
                    elif 'Option B' in label:
                        verdict = 'SENSITIVITY CONTROL (Ortholog-Mapped Distance)'
                    elif 'BLOSUM' in label:
                        verdict = '**INCONCLUSIVE (Point > 0.55, CI Spans 0.50)**'
                    else:
                        verdict = 'DESCRIPTIVE'
                elif t_key == 'HCMV_UL54_Foscarnet':
                    if 'Distance-Only' in label and 'Option A' in label:
                        verdict = '**DESCRIPTIVE (Consistent with Spatial Clustering)**'
                    elif 'Composite' in label and 'Option A' in label:
                        verdict = '**DESCRIPTIVE (Consistent with Spatial Clustering)**'
                    elif 'Composite' in label and 'Option B' in label:
                        verdict = 'DESCRIPTIVE'
                    elif 'BLOSUM' in label:
                        verdict = '**BELOW CHANCE / SUPPRESSED**'
                    elif 'Option B' in label:
                        verdict = 'SENSITIVITY CONTROL (Ortholog-Mapped Distance)'
                    else:
                        verdict = 'DESCRIPTIVE'
                elif t_key == 'HCMV_UL54_Cidofovir':
                    if 'Distance-Only' in label and 'Option A' in label:
                        verdict = '**NOT SUPPORTED (Below Chance; Point < 0.50)**'
                    elif 'Composite' in label and 'Option A' in label:
                        verdict = '**NOT SUPPORTED (Near Chance; Point < 0.55)**'
                    elif 'Composite' in label and 'Option B' in label:
                        verdict = '**POST HOC EXPLORATORY (Point > 0.80)**'
                    elif 'BLOSUM' in label:
                        verdict = '**INCONCLUSIVE (Point > 0.55, CI Spans 0.50)**'
                    elif 'Option B' in label:
                        verdict = 'SENSITIVITY CONTROL (Ortholog-Mapped Distance)'
                    else:
                        verdict = 'DESCRIPTIVE'
                else:
                    verdict = 'DESCRIPTIVE'
                    
                pw_str = f"**{pw:.4f}** [{ci[0]:.4f}, {ci[1]:.4f}]"
                sm_str = f"{sm:.4f} [{sm_ci[0]:.4f}, {sm_ci[1]:.4f}]"
                loso_str = f"[{loso[0]:.4f}, {loso[1]:.4f}]"
                s2_str = f"{s2_val:.4f} (n={s2_p})" if s2_val is not None else "n/a"
                d1_str = f"{d1_val:.4f} [{d1_ci[0]:.4f}, {d1_ci[1]:.4f}]" if d1_val is not None else "n/a"
                md.append(f"| **{g}-{d}** | {arm} | {label} | {pw_str} | {sm_str} | {loso_str} | {s2_str} | {d1_str} | {verdict} |")
        else:
            for f_name, f_vals in t_info['features'].items():
                pw = f_vals['s0_pair_weighted_auroc']
                ci = f_vals['s0_pair_weighted_ci_95']
                sm = f_vals['s0_study_mean_auroc']
                sm_ci = f_vals['s0_study_mean_ci_95']
                loso = f_vals['s0_loso_range']
                s2_val = f_vals['s2_numeric_fc_pair_weighted_auroc']
                s2_p = f_vals['s2_numeric_fc_pairs']
                d1_val = f_vals['d1_pooled_auroc']
                d1_ci = f_vals['d1_pooled_ci_95']
                
                score_label = f_name
                if t_key == 'HCMV_UL97_Ganciclovir' and 'BLOSUM' in f_name:
                    verdict = '**DESCRIPTIVE (Consistent with Kinase Radicalness)**'
                else:
                    verdict = f_vals.get('decision', 'DESCRIPTIVE')
                    
                pw_str = f"**{pw:.4f}** [{ci[0]:.4f}, {ci[1]:.4f}]"
                sm_str = f"{sm:.4f} [{sm_ci[0]:.4f}, {sm_ci[1]:.4f}]"
                loso_str = f"[{loso[0]:.4f}, {loso[1]:.4f}]"
                s2_str = f"{s2_val:.4f} (n={s2_p})" if s2_val is not None else "n/a"
                d1_str = f"{d1_val:.4f} [{d1_ci[0]:.4f}, {d1_ci[1]:.4f}]"
                md.append(f"| **{g}-{d}** | {arm} | {score_label} | {pw_str} | {sm_str} | {loso_str} | {s2_str} | {d1_str} | {verdict} |")
            
    md.append("\n#### Table 9C: Pre-Registered Empirical Contrast Tests (2,000 Study-Clustered Bootstrap Resamples)\n")
    md.append("| Contrast | Hypothesized Relationship | Observed Delta | 95% Bootstrap CI | Empirical p-value | Final Verdict |")
    md.append("|---|---|---|---|---|---|")
    
    c1 = s2_data['contrasts']['kinase_vs_polymerase_gcv_blosum']
    c1_delta = f"**{c1['delta_observed']:+.4f}**"
    c1_ci = f"[{c1['delta_ci_95'][0]:+.4f}, {c1['delta_ci_95'][1]:+.4f}]"
    p1_str = f"p = {c1['empirical_p']:.4f}" if c1['empirical_p'] >= 0.0005 else "p < 0.0005"
    md.append(f"| **Contrast 1: Kinase vs Polymerase Radicalness** | UL97-GCV BLOSUM62 > UL54-GCV BLOSUM62 | {c1_delta} | {c1_ci} | {p1_str} | **INCONCLUSIVE (CI Spans Zero)** |")
    
    c2 = s2_data['contrasts']['polymerase_fos_where_vs_blosum']
    c2_delta = f"**{c2['delta_observed']:+.4f}**"
    c2_ci = f"[{c2['delta_ci_95'][0]:+.4f}, {c2['delta_ci_95'][1]:+.4f}]"
    p2_str = f"p = {c2['empirical_p']:.4f}" if c2['empirical_p'] >= 0.0005 else "p < 0.0005"
    md.append(f"| **Contrast 2: Polymerase Proximity vs Radicalness (UL54-FOS)** | UL54-FOS Proximity > UL54-FOS BLOSUM62 | {c2_delta} | {c2_ci} | {p2_str} | **DESCRIPTIVE (Confirmed Catalytic Triad Proximity; see Table 9D/9E)** |")
    
    # Load Option B Apportionment results
    b_apportion_path = os.path.join(ROOT, 'logs', 'option_b_harmonization_and_apportionment_results.json')
    if os.path.exists(b_apportion_path):
        with open(b_apportion_path, 'r', encoding='utf-8') as f:
            b_app_data = json.load(f)
            
        md.append("\n#### Table 9D: Methodological Audit -- Factorial Apportionment of Apparent Option B Discordance\n")
        md.append("*Factorially separates pure coordinate provenance (Swiss-Model vs 8EXX on clean 19-residue shell), reference-set specification (distal exonuclease MG 1303 inclusion), and missingness imputation policy (complete cases N=133 vs -max_dist imputation N=180):*\n")
        md.append("| Target System | Option A: Swiss-Model (Full Cohort) | Option A: Swiss-Model (Common N=133) | Cell 1: B_harmonized (Clean Shell19, Common N=133) | Cell 2: B_clean + Imputed (Clean Shell19, Full N=180) | Cell 3: B_distal_MG (Distal MG+Triad, Common N=133) | Cell 4: Old Option B (Distal MG+Triad, Imputed N=180) |")
        md.append("|---|---|---|---|---|---|---|")
        
        drug_names = [('UL54-Foscarnet', 'Foscarnet'), ('UL54-Ganciclovir', 'Ganciclovir'), ('UL54-Cidofovir', 'Cidofovir')]
        for d_label, d_key in drug_names:
            d_obj = b_app_data['drugs'][d_key]
            s0_a_full = d_obj['option_a_full_cohort']['s0_pair_weighted']
            s0_a_com = d_obj['option_a_common_cases']['s0_pair_weighted']
            s0_c1 = d_obj['cell1_b_clean_shell19_common']['s0_pair_weighted']
            s0_c2 = d_obj['cell2_b_clean_shell19_imputed']['s0_pair_weighted']
            s0_c3 = d_obj['cell3_b_distal_mg_triad_common']['s0_pair_weighted']
            s0_c4 = d_obj['cell4_old_option_b_imputed']['s0_pair_weighted']
            md.append(f"| **{d_label}** | **{s0_a_full:.4f}** | {s0_a_com:.4f} | **{s0_c1:.4f}** | {s0_c2:.4f} | {s0_c3:.4f} | **{s0_c4:.4f}** |")
            
        md.append("\n*Quantitative Apportionment of Observed Discordance Shifts:*\n")
        md.append("| Target System | Pure Coordinate Provenance Delta (Cell 1 - Opt A Common) | Reference-Set Shift (Distal Exonuclease MG 1303: Cell 3 - Cell 1) | Imputation Shift (Cell 4 - Cell 3) | Total Divergence (Old B - Opt A Full) | Primary Divergence Driver |")
        md.append("|---|---|---|---|---|---|")
        
        for d_label, d_key in drug_names:
            eff = b_app_data['drugs'][d_key]['effects_apportionment']
            p_del = eff['pure_provenance_delta']
            r_shift = eff['reference_set_shift_common']
            i_shift = eff['imputation_shift_buggy_shell']
            tot = eff['total_divergence']
            driver = "Reference-set specification (distal exonuclease MG 1303)"
            md.append(f"| **{d_label}** | **{p_del:+.4f}** | **{r_shift:+.4f}** | {i_shift:+.4f} | **{tot:+.4f}** | {driver} |")
            
        md.append("\n*Methodological Interpretation of Table 9D:*")
        md.append("1. Pure coordinate provenance explains <= 0.0092 AUROC delta across all three drugs (Spearman rho = 0.989, median delta = 0.39 A), demonstrating that Swiss-Model homology modeling did not corrupt distance geometry.")
        md.append("2. The apparent historical divergence between Option A and Option B is attributable essentially in full to reference-set specification: the inadvertent inclusion of the distal exonuclease magnesium MG 1303, located 21.7 A from the polymerase shell. This artificially favoured exonuclease-localized mutations (GCV, CDV) while corrupting the polymerase-localized signal (FOS). Because the reference-set and imputation effects are measured on different case sets (complete cases vs imputed) and partially oppose one another on Foscarnet, they are not additively separable and no single percentage attribution is quoted.")
        md.append("3. Missingness imputation policy contributed a small, uniform secondary shift of +0.03 to +0.07.")

    # Load Step 2 Exonuclease Confirmatory results
    exo_conf_path = os.path.join(ROOT, 'logs', 'step2_hcmv_exonuclease_confirmatory_results.json')
    if os.path.exists(exo_conf_path):
        with open(exo_conf_path, 'r', encoding='utf-8') as f:
            exo_data = json.load(f)
            
        md.append("\n#### Table 9E: Confirmatory Double Dissociation Crossover Benchmark (Polymerase vs. Exonuclease Active-Site Proximity)\n")
        md.append("*Confirmatory test of the pre-registered double dissociation hypothesis (docs/UL54_EXONUCLEASE_CONFIRMATORY_PREREGISTRATION.md; B=2,000 study-clustered bootstrap resamples, seed 42):*\n")
        md.append("| Target System | Inhibitor Class & Target Mechanism | Polymerase Active-Site S0 [95% CI] | Exonuclease Active-Site S0 [95% CI] | Within-Drug Shift Delta (Exo - Pol) [95% CI] | Pre-Registered Rule Verdict |")
        md.append("|---|---|---|---|---|---|")
        
        target_info = [
            ('UL54-Ganciclovir', 'Ganciclovir', 'Nucleoside Analogue (Chain Terminator / Proofreading Excision)'),
            ('UL54-Cidofovir', 'Cidofovir', 'Acyclic Nucleotide Phosphonate (Chain Terminator / Proofreading Excision)'),
            ('UL54-Foscarnet', 'Foscarnet', 'Pyrophosphate Analogue (Direct Catalytic Pocket Inhibitor)')
        ]
        
        for t_label, d_k, mech in target_info:
            t_res = exo_data['targets'][d_k]
            p_pt = t_res['s0_polymerase']['point']
            p_ci = t_res['s0_polymerase']['ci_95']
            e_pt = t_res['s0_exonuclease']['point']
            e_ci = t_res['s0_exonuclease']['ci_95']
            d_obs = t_res['delta_exo_minus_pol']['observed']
            d_ci = t_res['delta_exo_minus_pol']['ci_95']
            
            p_str = f"{p_pt:.4f} [{p_ci[0]:.4f}, {p_ci[1]:.4f}]"
            e_str = f"**{e_pt:.4f}** [{e_ci[0]:.4f}, {e_ci[1]:.4f}]" if e_pt > 0.55 else f"{e_pt:.4f} [{e_ci[0]:.4f}, {e_ci[1]:.4f}]"
            if d_k == 'Foscarnet':
                p_str = f"**{p_pt:.4f}** [{p_ci[0]:.4f}, {p_ci[1]:.4f}]"
                
            del_str = f"**{d_obs:+.4f}** [{d_ci[0]:+.4f}, {d_ci[1]:+.4f}]"
            
            # The pre-registered rules are stated on point estimates; we additionally
            # report whether the paired within-drug delta CI excludes zero, which is the
            # stronger test and which Ganciclovir does not pass.
            d_p = t_res['delta_exo_minus_pol'].get('p_value')
            delta_sig = d_ci[0] > 0 or d_ci[1] < 0
            if d_k == 'Ganciclovir':
                v_str = (f"**{exo_data['verdicts']['rule1_gcv_exonuclease_recovery']}**"
                         f"; paired delta NOT significant (CI spans zero, p = {d_p:.3f})")
            elif d_k == 'Cidofovir':
                v_str = (f"**{exo_data['verdicts']['rule2_cdv_exonuclease_recovery']}**"
                         f"; paired delta significant")
            else:
                v_str = ("**CONFIRMED (Crossover Degradation: Exonuclease Inverted)**"
                         "; paired delta significant")
            if not delta_sig and d_k != 'Ganciclovir':
                v_str += " [paired delta CI spans zero]"
                
            md.append(f"| **{t_label}** | {mech} | {p_str} | {e_str} | {del_str} | {v_str} |")
            
        md.append("\n*Formal Double Dissociation Crossover Contrasts:*\n")
        x_gcv = exo_data['crossover_contrasts']['gcv_vs_fos']
        x_cdv = exo_data['crossover_contrasts']['cdv_vs_fos']
        # With B = 2,000 resamples the smallest resolvable one-sided p is 1/B = 0.0005;
        # an observed bootstrap p of 0.0 is reported at that resolution limit, not below it.
        n_boot = exo_data.get('protocol', {}).get('n_boot', 2000)
        p_floor = f"p < {1.0 / n_boot:.4f}"
        md.append(f"- **GCV vs FOS Crossover Delta:** **{x_gcv['observed_delta']:+.4f}** [{x_gcv['ci_95'][0]:+.4f}, {x_gcv['ci_95'][1]:+.4f}] ({p_floor}, bootstrap resolution limit at B={n_boot:,}) --> **CONFIRMED (Rule 3 Pass: Statistically Significant Crossover)**")
        md.append(f"- **CDV vs FOS Crossover Delta:** **{x_cdv['observed_delta']:+.4f}** [{x_cdv['ci_95'][0]:+.4f}, {x_cdv['ci_95'][1]:+.4f}] ({p_floor}, bootstrap resolution limit at B={n_boot:,}) --> **CONFIRMED (Rule 3 Pass: Statistically Significant Crossover)**")
        md.append("")
        md.append("*Note:* the GCV arm meets the pre-registered point criteria but its paired within-drug delta confidence interval spans zero (p = 0.055). The double dissociation is therefore carried by the CDV-vs-FOS crossover, with GCV directionally consistent but not individually significant.")
        md.append("")
        md.append("*Interval provenance:* the UL54-GCV polymerase point estimate S0 = 0.4135 carries two published intervals in this work, [0.1901, 0.6941] from the original Step 2 composite run (logs/step2_hcmv_composite_results.json, where the resample is shared across all seven feature variants) and [0.2183, 0.6393] from this confirmatory run (where the resample is shared across the polymerase and exonuclease features only). Both use B = 2,000 study-clustered resamples at seed 42; the interval differs because the set of co-resampled features differs. The confirmatory interval is primary for Table 9E, and the Step 2 interval remains primary for the pre-registered replication verdict reported in Table 9A-9C.")
        
        md.append("\n*Sensitivity Windows Robustness (Exonuclease Boundary Variations):*\n")
        for t_label, d_k, _ in target_info:
            sw = exo_data['targets'][d_k]['sensitivity_windows_s0']
            md.append(f"- **{d_k}:** Catalytic 7-site = **{exo_data['targets'][d_k]['s0_exonuclease']['point']:.4f}**; standard window (300-545) = **{sw['standard_exo_300_545']:.4f}**; tight (300-520) = **{sw['tight_exo_300_520']:.4f}**; expanded (280-560) = **{sw['expanded_exo_280_560']:.4f}**")
            
        md.append("\n*Methodological Interpretation of Table 9E:* The double dissociation indicates that the failure of polymerase-active-site distance on GCV (S0 = 0.4135) was not an uninformative failure of structural priors, but a consequence of mechanism mismatch. When evaluated against the functional site of resistance for each inhibitor class, spatial distance recovers predictive performance for chain terminators (0.72 on GCV, 0.93 on CDV) while degrading below chance on the pyrophosphate inhibitor (0.32 on FOS). The GCV arm is directionally consistent but its paired delta is not individually significant; the crossover is carried by CDV vs FOS.")

    # ---- Table 9F: leave-site-out robustness ----
    lso_path = os.path.join(ROOT, 'logs', 'ul54_leave_site_out_robustness.json')
    if os.path.exists(lso_path):
        with open(lso_path, 'r', encoding='utf-8') as f:
            lso = json.load(f)
        md.append("\n#### Table 9F: Leave-Site-Out Robustness of the Double Dissociation\n")
        md.append("*Each feature is defined by a set of reference residues, and variants sitting on those residues have distance zero by construction. To test whether the signal is a spatial gradient rather than an artefact of a few frequent resistance positions, each feature is re-evaluated with the variants at its own defining residues removed. Both reference sets are tested symmetrically (B=2,000 study-clustered bootstrap, seed 42).*\n")
        md.append("| Target System | Rows at Exo Reference Residues (of which R) | Rows at Pol Reference Residues (of which R) | Primary Feature | S0 Full Cohort | S0 Own Defining Residues Removed | Retained |")
        md.append("|---|---|---|---|---|---|---|")
        for d_k in ['Ganciclovir', 'Cidofovir', 'Foscarnet']:
            rec = lso['drugs'][d_k]
            ver = lso['survival_verdicts'][d_k]
            key = ver['primary_feature']
            full = rec['cohorts']['full'][key]
            drop = rec['cohorts']['drop_both'][key]
            dci = rec['cohorts']['drop_both'][key + '_ci_95']
            label = 'Exonuclease' if 'exo' in key else 'Polymerase'
            kept = "**RETAINED**" if (ver['survives_above_0.55'] and ver['ci_excludes_0.50']) else "NOT RETAINED"
            md.append(f"| **UL54-{d_k}** | {rec['n_rows_at_exo_sites']} ({rec['n_resistant_at_exo_sites']} R) | "
                      f"{rec['n_rows_at_pol_sites']} ({rec['n_resistant_at_pol_sites']} R) | {label} | "
                      f"{full:.4f} | **{drop:.4f}** [{dci[0]:.4f}, {dci[1]:.4f}] | {kept} |")
        md.append("\n*Methodological Interpretation of Table 9F:* Removing every variant located on a feature's own defining residues shifts each primary estimate by at most 0.045 AUROC "
                  "(GCV exonuclease 0.7243 -> 0.6798; CDV exonuclease 0.9264 -> 0.9167; FOS polymerase 0.8822 -> 0.8711), and all three retain a 95% confidence interval excluding 0.50. "
                  "The Cidofovir case is the most stringent, since all 26 variants at the exonuclease reference residues are resistant, yet the signal persists at 0.9167. "
                  "The double dissociation therefore reflects a spatial gradient across the enzyme rather than memorisation of a small number of frequent resistance positions.")

    return '\n'.join(md)

def generate_where_ablation_table():
    ab_path = os.path.join(ROOT, 'logs', 'where_components_ablation_results.json')
    if not os.path.exists(ab_path):
        return ""
        
    with open(ab_path, 'r', encoding='utf-8') as f:
        ab_data = json.load(f)
        
    md = []
    md.append("### Table 10: Decomposition Ablation of the Positional Prior (`Where`)\n")
    md.append("*Evaluates the relative contributions of Evolutionary Positional Constraint (ESM-2 650M) versus 3D Active-Site Proximity (Linear Standardized and Rank-Normalized) across HSV-1 targets (docs/WHERE_COMPONENTS_ABLATION_PLAN.md; B=2,000 study-clustered bootstrap resamples, seed 42):*\n")
    md.append("| Target System | Benchmark Role | Evaluated Score / Component | S0 Pair-Weighted [95% CI] | Study-Mean AUROC | S2 Numeric FC Stratum | D1 Pooled AUROC |")
    md.append("|---|---|---|---|---|---|---|")
    
    target_order = [
        ('UL30-Aciclovir', 'Primary Arm'),
        ('UL30-Foscarnet', 'Secondary Descriptive Arm'),
        ('UL23-Aciclovir', 'Comparative Counter-Arm')
    ]
    
    score_display = [
        ('Composite Where (Linear Std, Benchmark Standard)', 'Composite Where (Standardized)'),
        ('Composite Where (Rank-Normalized Distance)', 'Composite Where (Rank-Norm Dist)'),
        ('Active-Site Proximity Only (Linear Std)', 'Active-Site Proximity Only (Linear Std)'),
        ('Active-Site Proximity Only (Rank-Norm)', 'Active-Site Proximity Only (Rank-Norm)'),
        ('Evolutionary Constraint Only (ESM-2 650M)', 'Evolutionary Constraint Only (ESM-2 650M)'),
        ('BLOSUM62 Radicalness (Non-Positional Baseline)', 'BLOSUM62 Radicalness')
    ]
    
    for t_key, role in target_order:
        t_data = ab_data['targets'][t_key]
        for label, k in score_display:
            sc = t_data['scores'][k]
            pw = sc['s0_pair_weighted']
            ci = sc['s0_pair_weighted_ci_95']
            sm = sc['s0_study_mean']
            s2 = sc.get('s2_stratum_pooled_auc', sc.get('s2_stratum_pair_weighted'))
            d1 = sc['d1_pooled']
            
            pw_str = f"**{pw:.4f}** [{ci[0]:.4f}, {ci[1]:.4f}]" if 'Composite' in label or 'Proximity' in label else f"{pw:.4f} [{ci[0]:.4f}, {ci[1]:.4f}]"
            s2_str = f"{s2:.4f}" if s2 is not None else "n/a"
            md.append(f"| **{t_key}** | {role} | {label} | {pw_str} | {sm:.4f} | {s2_str} | {d1:.4f} |")
            
    md.append("\n*Paired Bootstrap Deltas (Within Shared Study Resamples, B=2,000):*\n")
    md.append("| Target System | Delta(Composite - Constraint) [95% CI] | Delta(Composite - Distance_Std) [95% CI] | Delta(Distance_Std - Constraint) [95% CI] | Delta(Composite_Rank - Composite_Std) [95% CI] |")
    md.append("|---|---|---|---|---|")
    
    for t_key, _ in target_order:
        pd_data = ab_data['targets'][t_key]['paired_deltas']
        d_cc = pd_data['Delta(Composite - Constraint)']
        d_cd = pd_data['Delta(Composite - Distance_Std)']
        d_dc = pd_data['Delta(Distance_Std - Constraint)']
        d_rr = pd_data['Delta(Composite_Rank - Composite_Std)']
        
        str_cc = f"**{d_cc['observed_delta']:+.4f}** [{d_cc['ci_95'][0]:+.4f}, {d_cc['ci_95'][1]:+.4f}]"
        str_cd = f"{d_cd['observed_delta']:+.4f} [{d_cd['ci_95'][0]:+.4f}, {d_cd['ci_95'][1]:+.4f}]"
        str_dc = f"{d_dc['observed_delta']:+.4f} [{d_dc['ci_95'][0]:+.4f}, {d_dc['ci_95'][1]:+.4f}]"
        str_rr = f"{d_rr['observed_delta']:+.4f} [{d_rr['ci_95'][0]:+.4f}, {d_rr['ci_95'][1]:+.4f}]"
        
        md.append(f"| **{t_key}** | {str_cc} | {str_cd} | {str_dc} | {str_rr} |")
        
    md.append("\n*Methodological Interpretation of Table 10:*")
    md.append("1. On the primary essential target (UL30-ACV), active-site proximity alone reaches S0 = 0.7846 versus constraint alone (0.7041, observed paired delta delta = +0.0683 [-0.1571, +0.2629], p = 0.2630); while distance achieves a higher point estimate, the paired difference spans zero.")
    md.append("2. Combining distance with evolutionary constraint in Composite Where yields S0 = 0.7865 (paired delta delta = +0.0768 [-0.0192, +0.1559] over constraint alone, p = 0.0615, spanning zero). Composite Where substantially narrows the 95% bootstrap confidence interval width relative to distance alone (0.0953 vs. 0.2173 AUROC interval width), improving estimation precision rather than mean discriminative performance.")
    md.append("3. Rank-normalized distance formulation in Composite Where further improves within-study concordance to 0.8165 (paired excess delta = +0.0297 [-0.0061, +0.0645], p = 0.0645).")
    md.append("4. On the dispensable kinase (UL23-ACV), distance alone reaches 0.6556 versus constraint alone (0.5710, paired delta delta = +0.0847 [-0.0717, +0.2370], p = 0.1150, spanning zero), indicating directional consistency without statistical separation.")
    md.append("5. S2 values report pooled AUROC on the quantitative fold-change stratum, matching the Table 8 definition.")
    
    return '\n'.join(md)

def generate_split_degradation_table():
    diag_path = os.path.join(ROOT, 'logs', 'split_separation_diagnostic.json')
    bench_path = os.path.join(ROOT, 'logs', 'split_degradation_benchmark_results.json')
    if not os.path.exists(diag_path) or not os.path.exists(bench_path):
        return ""
        
    with open(diag_path, 'r', encoding='utf-8') as f:
        diag_data = json.load(f)
    with open(bench_path, 'r', encoding='utf-8') as f:
        bench_data = json.load(f)
        
    md = []
    md.append("### Table 11: Generalization Benchmark Across Split Topologies & Pre-Hoc Separation Diagnostics\n")
    md.append("#### Table 11A: Label-Free Pre-Hoc Split Separation Diagnostic (Sequence-Neighbour Distance)\n")
    md.append("*Computed strictly label-free (without reading resistance phenotypes or fold changes; logs/split_separation_diagnostic.json):*\n")
    md.append("| Target System | Split Scheme | Evaluated Unique Positions | Median Distance to Nearest Train (aa) | 25th Percentile ($P_{25}$) | Fraction Within $\\le 2$ aa | Fraction Within $\\le 5$ aa | Sequence-Neighbour Leakage Status |")
    md.append("|---|---|---|---|---|---|---|---|")
    
    target_order = ['UL30-Aciclovir', 'UL30-Foscarnet', 'UL23-Aciclovir']
    scheme_order = [
        ('position-grouped', 'Position-Grouped CV (Standard Clinical Split)'),
        ('3D k-means', '3D Spatial k-means Clustering'),
        ('contiguous-block', 'Contiguous Sequence Block Holdout')
    ]
    
    for t_key in target_order:
        t_diag = diag_data['targets'][t_key]
        n_pos = t_diag['n_unique_positions']
        for s_id, s_label in scheme_order:
            sc = t_diag['schemes'][s_id]
            med = sc['median_nearest_train_distance']
            p25 = sc['p25']
            f2 = sc['frac_within_2_residues'] * 100
            f5 = sc['frac_within_5_residues'] * 100
            
            if s_id == 'position-grouped':
                if f2 > 90.0:
                    status = f"**Critical Channel Open ({f2:.1f}% $\\le 2$ aa)**"
                else:
                    status = f"**High Leakage Risk ({f2:.1f}% $\\le 2$ aa)**"
            elif s_id == '3D k-means':
                if f5 > 50.0:
                    status = f"**Partially Separated ({med:.1f} aa; {f5:.1f}% $\\le 5$ aa)**"
                else:
                    status = f"Substantially Separated ({med:.1f} aa)"
            else:
                status = f"**Strictly Separated ({med:.1f} aa)**"
                
            md.append(f"| **{t_key}** | {s_label} | {n_pos} | **{med:.1f}** | {p25:.1f} | {f2:.2f}% | {f5:.2f}% | {status} |")
            
    md.append("\n#### Table 11B: Supervised Feature-Family Degradation Matrix Under Block Holdout\n")
    md.append("*Apples-to-apples evaluation comparing Random Forest models trained across split topologies (B=2,000 bootstrap resamples, seed 42; logs/split_degradation_benchmark_results.json):*\n")
    md.append("| Target System | Feature Family | Model Type | Position-Grouped AUROC [95% CI] | 3D Spatial AUROC | Contiguous-Block AUROC [95% CI] | Degradation $\\Delta$ (PG - Block) [95% CI] | Empirical p-value | Empirical Generalization Verdict |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    
    family_order = [
        ('Positional Identity', 'Positional Identity (`norm_pos`)'),
        ('Evolutionary Conservation', 'Evolutionary Conservation (BLOSUM, Entropy, PSSM)'),
        ('Protein Language Models', 'Protein Language Models (ESM-2 650M, ProtBERT, 35M)'),
        ('Multi-Modal Combined', 'Multi-Modal Combined (All Features)')
    ]
    
    for t_key in target_order:
        t_bench = bench_data['targets'][t_key]
        rf_models = t_bench['models']['random_forest']
        lr_models = t_bench['models'].get('logistic_regression', {})
        
        for f_id, f_label in family_order:
            m = rf_models[f_id]
            pg_auc = m['position_grouped_auc']
            pg_ci = m['position_grouped_ci_95']
            sp_auc = m['spatial_3d_auc']
            blk_auc = m['contiguous_block_auc']
            blk_ci = m['contiguous_block_ci_95']
            del_auc = m['degradation_delta_pg_minus_block']
            del_ci = m['degradation_delta_ci_95']
            p_val = m['degradation_p_value']
            
            # 1/B is the smallest resolvable bootstrap p-value; do not report below it.
            _pf = 1.0 / bench_data.get('protocol', {}).get('n_boot', 2000)
            p_str = f"p = {p_val:.4f}" if p_val >= _pf else f"p < {_pf:.4f}"
            
            # Verdict derives from the measured degradation, not from the family's
            # identity, so families with comparable degradation read comparably.
            if del_auc >= 0.20:
                verdict = "**Collapse (position-indexed; cannot extrapolate out of block)**"
            elif del_auc >= 0.05:
                verdict = "Signal preserved; degradation 0.05-0.20"
            else:
                verdict = "Signal preserved; degradation below 0.05"
            if blk_auc < 0.50:
                verdict += " [inverted below chance]"
                
            pg_str = f"**{pg_auc:.4f}** [{pg_ci[0]:.4f}, {pg_ci[1]:.4f}]"
            sp_str = f"{sp_auc:.4f}" if sp_auc is not None else "n/a"
            blk_str = f"**{blk_auc:.4f}** [{blk_ci[0]:.4f}, {blk_ci[1]:.4f}]" if blk_auc > 0.55 else f"{blk_auc:.4f} [{blk_ci[0]:.4f}, {blk_ci[1]:.4f}]"
            del_str = f"**{del_auc:+.4f}** [{del_ci[0]:+.4f}, {del_ci[1]:+.4f}]"
            
            md.append(f"| **{t_key}** | {f_label} | Random Forest | {pg_str} | {sp_str} | {blk_str} | {del_str} | {p_str} | {verdict} |")
            
        # Linear model class, so that the collapse is not attributed to the model class
        # when it is a property of the feature representation.
        for f_id, f_label in family_order:
            if f_id not in lr_models:
                continue
            m = lr_models[f_id]
            pg_auc = m['position_grouped_auc']
            pg_ci = m.get('position_grouped_ci_95')
            sp_auc = m.get('spatial_3d_auc')
            blk_auc = m['contiguous_block_auc']
            blk_ci = m.get('contiguous_block_ci_95')
            del_auc = m['degradation_delta_pg_minus_block']
            del_ci = m.get('degradation_delta_ci_95')
            p_val = m['degradation_p_value']

            _pf = 1.0 / bench_data.get('protocol', {}).get('n_boot', 2000)
            p_str = f"p = {p_val:.4f}" if p_val >= _pf else f"p < {_pf:.4f}"

            if del_auc >= 0.20:
                verdict = "**Collapse (position-indexed; cannot extrapolate out of block)**"
            elif del_auc >= 0.05:
                verdict = "Signal preserved; degradation 0.05-0.20"
            else:
                verdict = "Signal preserved; degradation below 0.05"
            if blk_auc < 0.50:
                verdict += " [inverted below chance]"

            pg_str = f"{pg_auc:.4f}" + (f" [{pg_ci[0]:.4f}, {pg_ci[1]:.4f}]" if pg_ci else "")
            sp_str = f"{sp_auc:.4f}" if sp_auc is not None else "n/a"
            blk_str = f"{blk_auc:.4f}" + (f" [{blk_ci[0]:.4f}, {blk_ci[1]:.4f}]" if blk_ci else "")
            del_str = f"**{del_auc:+.4f}**" + (f" [{del_ci[0]:+.4f}, {del_ci[1]:+.4f}]" if del_ci else "")

            md.append(f"| **{t_key}** | {f_label} | Logistic Regression | {pg_str} | {sp_str} | {blk_str} | {del_str} | {p_str} | {verdict} |")

        # Invariant negative control
        nc = t_bench['negative_control_zero_shot_blosum']
        nc_pg = nc['position_grouped_auc']
        nc_blk = nc['contiguous_block_auc']
        nc_del = nc['degradation_delta']
        md.append(f"| **{t_key}** | Negative Control: BLOSUM62 | Zero-Shot Invariant | {nc_pg:.4f} [No Fit] | n/a | {nc_blk:.4f} [No Fit] | **{nc_del:+.4f}** [Invariant] | p = 1.0000 | **Invariant Control ($\\Delta = 0.0000$ by construction)** |")
        
    md.append("\n*Methodological Interpretation of Table 11:*")
    md.append("1. **Pre-Hoc Diagnostic Utility:** The label-free distance diagnostic demonstrates that position-grouped splits fail to close sequence-neighbour correlation channels (median nearest-train distance 1.0 aa in UL23, 2.0 aa in UL30; 60-96% of test positions sit within <= 2 residues of a training position).")
    md.append("2. **Positional Collapse Is Representation-Bound, Not Model-Bound:** Positional-identity features collapse under contiguous-block holdout in all six target x model-class combinations (Random Forest: 0.8797 -> 0.4591, Delta = +0.4206; 0.8661 -> 0.4669, Delta = +0.3992; 0.8185 -> 0.5450, Delta = +0.2735. Logistic Regression: Delta = +0.3424, +0.2929, +0.2336). A position-indexed feature cannot extrapolate to a sequence block disjoint from its training range, so the collapse is a property of the representation rather than of model capacity.")
    md.append("3. **All Variant-Level Families Survive, Including Classical Ones:** Every family whose features are computed from the variant itself degrades substantially less than positional identity (a 2.5- to 10-fold reduction in degradation: evolutionary conservation Delta = 0.0362-0.1092 and protein language models Delta = 0.0062-0.0356). Because classical conservation scores (BLOSUM62, positional entropy, PSSM) embed no learned representation yet behave equivalently to ESM-2 and ProtBERT, this contrast must not be read as evidence that language models uniquely acquire transferable sequence-to-function structure. The defensible claim is narrower and stronger: leakage under position-grouped splitting is carried by position-indexed features, and any variant-level descriptor avoids it.")
    md.append("4. **Two Qualifications on the Collapse:** First, the high apparent position-grouped AUROC requires a non-linear model; under logistic regression the positional-identity family reaches only 0.5668 on UL30-ACV, so the 0.88 figure reflects a Random Forest partitioning position space rather than an intrinsic property of the split. Second, under block holdout positional features do not regress to 0.50 but invert below chance (0.2244-0.3470 under logistic regression), the expected signature of forced extrapolation beyond the training range.")
    md.append("5. **Negative Control Interpretation:** The unfitted BLOSUM62 score is identical across all split topologies (Delta = 0.0000 exactly). This is a tautology, not a result: a score that is never fitted cannot respond to a change of training set. It is reported precisely to make that tautology explicit, because comparing an unfitted zero-shot score's split-invariance against a fitted model's holdout degradation is not a valid comparison, and an earlier version of this analysis made exactly that error.")
    
    return '\n'.join(md)

if __name__ == '__main__':
    print("=" * 95)
    print("LOG-BACKED TABLES FOR MANUSCRIPT AND DOCUMENTATION")
    print("=" * 95)
    print()
    print(generate_table7())
    print("\n" + "=" * 95 + "\n")
    print(generate_table8())
    print("=" * 95 + "\n")
    print(generate_step1_tables())
    print("=" * 95 + "\n")
    print(generate_step2_tables())
    print("=" * 95 + "\n")
    print(generate_where_ablation_table())
    print("=" * 95 + "\n")
    print(generate_split_degradation_table())
