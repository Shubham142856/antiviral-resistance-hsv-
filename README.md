# Position-only prediction in HSV-1 resistance records and a BRCA1 functional assay

The [research manuscript](manuscript/manuscript_draft.md) asks how well residue number alone distinguishes literature-curated phenotypes when a benchmark holds out individual positions, structural clusters or contiguous sequence regions.

On the same 262 structure-mapped UL30–aciclovir records, pooled AUROC was **0.875** for position grouping, **0.675** for structural clusters and **0.523** for sequence blocks. Three of four evaluable blocks had AUROC 0.500. The fifth held a single phenotype. On all 332 records in that arm, a control assigning one training-prevalence score per held-out block yielded pooled AUROC 0.190. Within every evaluable block it scored exactly 0.500. The manuscript decomposes the pooled comparisons to explain the difference.

## Reproduce the main result

Run the following command from the repository root. The verified CPU run used Python 3.13.14, NumPy 2.5.0, pandas 2.3.3, scikit-learn 1.9.0, SciPy 1.16.1 and Matplotlib 3.10.5. It reads the stored fold assignments.

```sh
python scripts/audit_position_control.py
```

The command generates the [fold-level results and paired intervals](logs/position_control_folds.json) and [manuscript figure](manuscript/figures/position_control_folds.png). The JSON includes the SHA-256 of the [input table](data/processed/primary_3d_spatial_folds.tsv) and the scikit-learn version. Each fold records the number of records, positive labels, unique positions and within-fold AUROC. A null AUROC identifies a single-class test fold.

The full dataset contains 929 HSV-1 missense–drug records across three arms. The three-way comparison uses the 737 records with resolved coordinates. The two-way sequence comparison also includes the remaining 192. Split construction is documented in the [position-fold builder](scripts/build_reproducible_folds_and_subsets.py) and [structural-fold builder](scripts/run_spatial_3d_holdout_experiment.py).

## Independent experimental assay

In 1,917 BRCA1 missense SNVs spanning 326 measured positions, pooled position-only AUROC fell from 0.674 to 0.449 when entire groups of neighbouring sites were held out. The paired difference was 0.225 (position-bootstrap 95% interval 0.135 to 0.311). Within an exon, adjacent positions differed by an average of 25.5 percentage points in their fraction of nonfunctional variants, compared with 29.5 points after shuffling sites within each exon. The endpoint is the authors' cellular assay classification.

The external audit uses two records of the Findlay et al. BRCA1 experiment. The [MaveDB CSV](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000097-0-2/scores) supplies normalized scores. The [authors' GEO supplement](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117159/suppl/GSE117159_Supplementary_Table_1.xlsx) supplies protein positions, exon experiments and functional classes. The script joins both exports by transcript variant and checks score agreement. The [analysis plan](docs/BRCA1_EXTERNAL_AUDIT.md) specifies the comparisons and decision criterion.

From this workspace in PowerShell, the following runs the external audit on CPU with temporary files and the plotting cache under the project directory.

```powershell
$auditCache = Join-Path (Get-Location) '.audit-cache'
New-Item -ItemType Directory -Force -Path $auditCache | Out-Null
$env:TMP = $auditCache
$env:TEMP = $auditCache
$env:MPLCONFIGDIR = $auditCache
python -B scripts/audit_brca1_external.py
```

The [external audit results](logs/brca1_external_audit.json) include data and protocol hashes, all fold denominators, the paired interval and the within-exon permutation test.

## Corrections to the earlier analysis

- The primary three-way comparison now uses identical records for every split. The primary arm has a median contiguous-block gap of 78.5 residues among coordinate-mapped records and 85 residues among all records.
- The historical six-target audit measured fold geometry. Under modulo-five assignment its observed range within two residues was 20.0% to 97.6%. [Geometry results](logs/cross_pathogen_leakage_audit.json)
- The historical clinical analysis mixed HSV types and used default feature values for missense records in its full-cohort calculation. The corrected accounting reuses the exact same missense classifications after verifying one-to-one record and phenotype agreement. Of 352 labelled HSV-1 UL23–aciclovir records, 351 qualify for this retrospective calculation and one combination record is excluded. PPV is 71.8% on those 351 records and 59.5% on the 286 missense records. The former 94.4% estimate is withdrawn. Clinical use would require independent validation. [Corrected counts](logs/clinical_triage_evaluation.json) · [Reproduction script](scripts/evaluate_clinical_triage_decision_rule.py)

Local Git revision `c14c30a` retains the earlier draft and analysis history. The repository also contains previously reported exploratory analyses. The current manuscript identifies the results used for this paper.

## Submission status

The project records three outstanding author-side steps: [dual review of the phenotype sample](logs/kappa_sample_50rows.tsv), [threshold citation verification](docs/threshold_citation_checklist.md) and [PI approval](docs/pi_signoff_form.md). Author names, affiliations, correspondence, funding and competing-interest declarations require confirmation before submission.
