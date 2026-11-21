# Findings from the full HSV antiviral resistance project

## Research question

How much evidence do published HSV-1 resistance records provide that sequence-based models generalize to previously uncharacterized regions, across publications, and to later reports? These tests address different aspects of generalization. Each requires its own denominator and comparison.

## Study population

The primary dataset comprises 929 single-site missense–drug records representing 630 distinct HSV-1 substitutions. UL30–aciclovir contributes 332 records, UL30–foscarnet 311 and UL23–aciclovir 286. The two UL30 arms share a protein and some substitutions. Literature curation determines which positions, drugs and phenotypes enter the dataset. The unit of analysis is a curated record, and the number of records cannot be interpreted as the number of independent patients. [Primary data](../data/processed/primary_3d_spatial_folds.tsv)

## Results

### Residue number establishes how much the partition matters

Across the same 262 coordinate-mapped UL30–aciclovir records, the position-only control scored 0.875 under position grouping, 0.675 under structural clustering and 0.523 under contiguous sequence holdout. The grouped-minus-block difference was 0.352 (conditional position-bootstrap 95% interval 0.253 to 0.447). On all 332 records, the pooled block AUROC was 0.459 although every evaluable block scored at least 0.500. Most positive–negative comparisons crossed fold boundaries. This control establishes the evaluation context for interpreting the model comparisons below. [Fold-level audit](../logs/position_control_folds.json)

### Model comparisons depend on the evaluation question

In the original within-block comparison of the supervised predictor against zero-shot ESM-2, the paired mean AUROC differences were +0.079 for UL23–aciclovir, −0.005 for UL30–aciclovir and +0.004 for UL30–foscarnet. All three block-level intervals included zero. The small number of evaluable blocks limits the precision of this comparison. The supervised predictor also achieved Brier scores of 0.178, 0.114 and 0.084 respectively on the evaluable-block cohorts. Each improved on a constant score set to that cohort's observed prevalence. The Brier score combines calibration and discrimination and therefore cannot establish calibration by itself. [Paired comparison and Brier audit](../logs/phase1_deliverables_results.json)

The exploratory feature-family comparison asked a different question. Both the protein-language-model features and the classical conservation features were supplied to fitted classifiers and compared on identical out-of-block records. On UL30–aciclovir, the pooled AUROC difference favouring the language-model family was 0.150 with a random forest (position-bootstrap 95% interval 0.070 to 0.232) and 0.215 with logistic regression (0.149 to 0.290). UL23–aciclovir gave inconclusive contrasts with both classifiers. The hypothesis was selected after inspecting these data. The two classifiers are sensitivity checks on one cohort. [Paired feature-family results](../logs/paired_family_contrast_results.json)

### Within-study comparisons retain a role for classical conservation

The evidence ladder evaluates pooled records, buffered regional holdouts, within-publication comparisons and a quantitative phenotype stratum. In the within-publication analysis, the classical evolutionary profile had pair-weighted concordance between 0.699 and 0.774 across the three arms. Every paired interval comparing that profile with ESM-2 included zero. This supports retaining the classical profile as a benchmark and leaves equivalence unresolved. For the primary arm, only 14 publications supplied both phenotypes, yielding 267 pairs. The paired advantage of the structural/evolutionary composite over the substitution baseline was 0.230 (study-clustered interval 0.015 to 0.496). The between-gene difference in that advantage was 0.141 (−0.118 to 0.482). Neither a cross-gene mechanism nor the individual contribution of each composite component is established by those estimates. [Evidence ladder](../logs/evidence_ladder_table8_results.json) · [Study-level controls](../logs/step1_ascertainment_control_results.json) · [Component ablation](../logs/where_components_ablation_results.json)

### A later publication can still describe a previously studied position

In the historical time split, the supervised AUROC on 131 later UL23–aciclovir records was 0.764. Among the 33 records at previously observed positions it was 0.887, versus 0.700 among 98 records at new positions. Later UL30 reports also included both seen and unseen positions. Each UL30 seen-position stratum contained only two positive records, making its AUROC particularly unstable. A publication-date split and a position holdout therefore answer distinct questions. These are retrospective subgroup estimates, with no prospective clinical endpoint. [Temporal confounding audit](../logs/temporal_confounding.log)

### Extensions and limitations

The HCMV comparison belongs in the account of completed work. The original cross-herpesvirus analysis encountered a replication failure, followed by a drug-dependent structural association on reanalysis. The hypothesis for the follow-up arose from the same HCMV data. Separate reporting of the initial result and the exploratory follow-up preserves that chronology. The cross-pathogen geometry audit likewise describes separation between folds. It did not test predictive performance in those other pathogens. [HCMV replication](../logs/step2_hcmv_replication_results.json) · [Cross-pathogen geometry audit](../logs/cross_pathogen_leakage_audit.json)

The historical clinical-rule analysis requires the corrected HSV-1 denominator. Among 352 labelled UL23–aciclovir records, 351 were evaluable under the rule and one combination record was excluded. The retrospective positive predictive value was 71.8% on the evaluable cohort and 59.5% on 286 missense records. The earlier 94.4% value mixed HSV types and used inconsistent inputs. The corrected estimates remain dependent on the curated cohort's case mix. [Corrected accounting](../logs/clinical_triage_evaluation.json)

## Interpretation

The original project's contribution is a comparison of what several evaluation designs allow researchers to conclude from the same sparse HSV literature. The position-only control quantifies regional interpolation. Paired model comparisons show where the choice of representation mattered in this dataset. The within-study and temporal analyses show how publication and site composition change the comparisons a reader can make. The full-project paper can report those findings together while keeping the BRCA1 external assay in its separate methods manuscript.

For supporting material, the original manuscript at Git revision `c14c30a` contains the zero-shot model tables, drug-conditioning and negative controls, structural coverage analysis, threshold and radius checks, and software description. Coverage gaps identify missing evidence. They do not establish why variants are absent. The open-source split diagnostic is a reproducibility contribution. [Sensitivity checks](../logs/phase3_verification_gates_results.json) · [Split diagnostic](../splitgap)
