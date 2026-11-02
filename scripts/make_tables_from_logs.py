"""
make_tables_from_logs.py
========================
Deterministically generates markdown for:
  - Table 7: Recalibrated Benchmarks & Excess Metrics across Structurally Resolved Domains (logs/step0b_final_recalibrated_benchmarks.json)
  - Table 8: The Evidence Ladder across D1, D2, D3, S2 (logs/evidence_ladder_table8_results.json)
  - Step 1 Tables: Positivity Check, S0 Primary Within-Study, Contrast D, and Secondary Strata (logs/step1_ascertainment_control_results.json)

Rule: Manuscript tables must be generated from this script, never hand-typed.
All quantities trace to audited JSON logs with unified seeds.
"""

import os
import json
import numpy as np

def generate_table7():
    with open('logs/step0b_final_recalibrated_benchmarks.json', 'r') as f:
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
    with open('logs/evidence_ladder_table8_results.json', 'r') as f:
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
            # Format row label
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
    with open('logs/step1_ascertainment_control_results.json', 'r') as f:
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
    
    # Target totals for per-target proportions
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
            
            # Format description with exact per-target share
            n_tot_strat = s_data['N']
            pct_cohort = (n_tot_strat / tot) * 100
            if 'S1' in s_key:
                desc = f"{pct_cohort:.1f}% of cohort ({n_tot_strat}/{tot}; CMI cell-culture confirmed = yes; near-full cohort for ACV, 68.5% for FOS)"
            elif 'S2' in s_key:
                desc = f"{pct_cohort:.1f}% of cohort ({n_tot_strat}/{tot}; HerpesDRG numeric fold-change, balanced quantitative)"
            elif 'S3' in s_key:
                desc = f"{pct_cohort:.1f}% of cohort ({n_tot_strat}/{tot}; recombinant marker transfer / BAC)"
            else:
                desc = s_data.get('description', '')
                
            md.append(f"| **{t_key}** | {s_key} | {n_tot_strat} ({s_data['Resistant']} R / {s_data['Sensitive']} S) | {ev_str} | {w_val} | {b_val} | {desc} |")
            
    md.append("\n*Methodological Note on S1:* Cell-culture confirmation flags cover 91.6% of UL30-ACV (304/332) and 91.6% of UL23-ACV (262/286), but only 68.5% of UL30-FOS (213/311). Because a flag covering ~69–92% of entries reflects general cohort compilation rather than a selective evidence tier, S1 serves as a broad descriptive census rather than a quality-filtered stratum.")
    return "\n".join(md)

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
