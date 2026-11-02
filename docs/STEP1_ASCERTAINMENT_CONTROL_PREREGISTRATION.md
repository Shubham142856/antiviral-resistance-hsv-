# Pre-Registration: Step 1 Ascertainment Bias & Curation Controls

**Document Version:** 1.0 (Frozen Pre-Registration)  
**Date:** 2026-11-02  
**Target Study:** HSV-1 Antiviral Resistance Benchmark & Mechanism Evaluation  
**Status:** Pre-registered BEFORE Execution. The Git commit SHA and timestamp serve as cryptographic proof of prior registration.  
**Repository:** [`https://github.com/Shubham142856/antiviral-resistance-hsv-`](https://github.com/Shubham142856/antiviral-resistance-hsv-)

---

## 1. Scientific Background & Rationale

A central vulnerability in retrospective viral resistance benchmarks is **curation ascertainment bias**:
If "sensitive" variants in public databases are predominantly distal, unmeasured natural polymorphisms cataloged during routine clinical sequencing, while "resistant" variants are proximal mutations characterized in targeted drug-pocket phenotyping studies, then any spatial or evolutionary positional prior ("Where") might achieve high apparent AUROC simply by separating surveillance polymorphisms from targeted drug-site studies.

To formally test whether the predictive signal of the zero-parameter `Where` prior reflects genuine biophysical fitness intolerance versus database curation artifacts, we evaluate the frozen `Where` metric within rigorous evidence strata and on quantitative phenotypic assay subsets.

---

## 2. Dataset Columns & Ground-Truth Definitions

Data is merged from the authoritative HerpesDRG database (`herpesdrg-db/herpesdrg-db.tsv`, commit `c7e2d9b`):

### Real HerpesDRG Header Columns:
```tsv
mutation_id	virus	genotype	gene	aa_change	Ganciclovir	Aciclovir	Cidofovir	Foscarnet	Brincidofovir	Letermovir	Brivudine	Penciclovir	Tomeglovir	Maribavir	Cyclopropavir	Amenamevir	Pritelivir	ref_title	ref_link	ref_doi	strain_generation	test_method	co_gene	co_aa	created_date	created_by	note	control_strain	status
```

### Verified Column Values for HSV-1 (`virus == 'HSV1'`):
1. **`strain_generation` (Evidence Tier):**
   - **Tier 1 (Recombinant / Lab-Engineered):** `Marker Transfer`, `Marker transfer`, `Recombinant BAC`, `Lab derived mutant`, `Lab induced strain`. These variants are engineered into defined viral reference backbones via marker transfer or site-directed mutagenesis, isolating the causal effect of the substitution independent of patient background.
   - **Tier 2 (Clinical Isolates):** `Isolated strain`. Naturally occurring patient isolates with verified susceptibility testing.
2. **Quantitative Phenotypic Assay Subset (`has_measured_fold_change`):**
   - The drug column (`Aciclovir` or `Foscarnet`) contains a direct **quantitative numeric IC50/EC50 fold-change** (e.g. `2.5`, `3.1`, `0.8`), verifying that the variant underwent quantitative laboratory testing.
   - Variants cataloged merely with categorical text (`Polymorphism` or `Sensitive`) without an experimental fold-change measurement are excluded from this stratum. Both resistant ($y=1$) and sensitive ($y=0$) arms must have measured fold-change data.
3. **`test_method` (Experimental Assay Standardization):**
   - **Plaque Reduction Assay (PRA):** Gold-standard phenotypic assay (`test_method == 'PRA'`).
   - **Dye Uptake Assay:** Standardized dye uptake spectrophotometric assay (`test_method == 'dye uptake assay'`).

---

## 3. Evaluated Model & Frozen Specifications

- **Evaluated Metric:** The frozen, zero-parameter spec-exact `Where` composite:
  $$\text{Where}(p) = Z(\text{Constraint}_{\text{PLM}}(p)) - Z(d_{\text{site}}(p))$$
  where $\text{Constraint}_{\text{PLM}} = -\mu_{\text{raw}}$ from ESM-2 (650M) 19-variant masked marginals, and $d_{\text{site}}$ is the Euclidean distance to the drug-binding pocket in `8EXX`/`8V1Q` (UL30) or `1OF1` (UL23).
- **No Refitting:** `Where` has zero learned weights. It is evaluated as a frozen scoring function.
- **Factor Decomposition:** We also evaluate the two individual components:
  1. Evolutionary Constraint alone: $Z(\text{Constraint}_{\text{PLM}})$
  2. Active-Site Proximity alone: $-Z(d_{\text{site}})$

---

## 4. Evaluated Target Arms

1. **UL30-Aciclovir ($N=332$, 62 resistant):** **Primary Hypothesis Target**.
2. **UL30-Foscarnet ($N=311$, 33 resistant):** Secondary Confirmatory Target.
3. **UL23-Aciclovir ($N=286$, 120 resistant):** Secondary Comparative Target.

---

## 5. Pre-Registered Decision Thresholds & Falsification Rules

For each evaluable stratum with $N_{\text{pos}} \ge 10$ events:
Let $\text{AUROC}_{\text{full}}$ be the full-cohort benchmark AUROC (UL30-ACV = 0.8606).
Let $\text{AUROC}_{\text{stratum}}$ be the within-stratum AUROC.

- **PASS (Genuine Biological Prior):**
  $$\text{AUROC}_{\text{stratum}} \ge \text{AUROC}_{\text{full}} - 0.05$$
  The positional prior maintains high discriminative power within the gold-standard recombinant and quantitatively measured subsets, confirming that ascertainment bias does not explain the signal.
- **INCONCLUSIVE:**
  $$\text{AUROC}_{\text{full}} - 0.10 \le \text{AUROC}_{\text{stratum}} < \text{AUROC}_{\text{full}} - 0.05$$
- **ARTIFACT (Curation Bias Driver):**
  $$\text{AUROC}_{\text{stratum}} < \text{AUROC}_{\text{full}} - 0.10$$
  If the within-stratum AUROC drops by $> 0.10$ upon excluding unmeasured polymorphisms or restricting to recombinant assays, the claim that `Where` is a true biological prior is **falsified**, and the signal will be explicitly reported as driven by curation ascertainment bias.

### Minimum Event Power Rule:
Any stratum with $N_{\text{pos}} < 10$ resistant events is acknowledged upfront as underpowered, reported descriptively, and skipped for binary pass/fail hypothesis decisions.
