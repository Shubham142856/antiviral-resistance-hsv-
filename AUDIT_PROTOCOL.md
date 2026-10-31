# AUDIT_PROTOCOL.md (complete text, version 4: strict)

**Project:** HSV antiviral-resistance mutation prediction (Paper 2)  
**Funding:** MPCST grant on computational approaches to HSV protein sequence analysis  
**Purpose:** Try to falsify each claim against the data, the power and the literature before any modeling. Every number marked "proposed" is mine and must be justified or replaced when you freeze the file.

**Go/no-go:** modeling starts only after `AUDIT_REPORT.md` exists, with a locked scope statement signed by the PI and a second annotator.

---

## 1. Stage 0: Frozen definitions

Commit this section, with its hash, before any counting. Ideally timestamp it on OSF or similar.

### 1.1 Naming
- ACV = acyclovir. **FOS = foscarnet** (never FCV, which is famciclovir, a TK-dependent prodrug). PTV = pritelivir. AMNV = amenamevir. CDV and others are not modeled.

### 1.2 Scope

| Tier | Gene | Drugs |
| :--- | :--- | :--- |
| Primary candidate | UL23 | ACV, FOS (FOS is a negative-control pair) |
| Primary candidate | UL30 | ACV, FOS |
| Provisional | UL5, UL52 | PTV, AMNV (primary only if the audit counts support it) |
| Secondary | UL8, UL29 | PTV/AMNV, ACV |
| Excluded | UL9, UL42, UL15, UL28 | no resistance evidence |

Record HSV-1 or HSV-2 on every row. Proposed: per-drug claims for HSV-1 only, with HSV-2 as a cross-type generalization test. PI to confirm.

### 1.3 Two separate fields (replaces the old A/B/C/D tiers)

| Field | Values | Meaning |
| :--- | :--- | :--- |
| `evidence_tier` | A: recombinant, BAC or marker-transfer isolating the single mutation. B: isolate or lab-derived virus where other changes may be present. C1: susceptibility individually phenotyped. C2: assumed from susceptible isolates, never individually phenotyped. D: no usable phenotype | How the phenotype was obtained |
| `phenotype_label` | Resistant, Susceptible, Polymorphism, Unknown | The outcome, assigned by the frozen threshold |

Separating them lets a recombinant mutation that turns out susceptible still have a tier.

### 1.4 Resistance thresholds
Fix per drug and assay a fold-change cutoff **with a cited source** (for example the CMI 2025 review, the HerpesDRG paper, CHARMD conventions, or an absolute IC50 cutoff). Rows from sources reporting only categories are recorded as reported and flagged. Threshold and source are stored per row. Placeholder table:

| Drug | Assay | Cutoff | Citation |
| :--- | :--- | :--- | :--- |
| ACV, FOS, PTV, AMNV | each assay type | TO BE CITED | TO BE CITED |

### 1.5 Temporal rule (choose one now)
- **Option A (natural split):** train on mutations in pre-2016 datasets; test on mutations newly reported in 2016 to 2024. The 2025 CMI review used a search window of January 2016 to September 2024, plus a pre-2016 dataset, and identified 221 novel resistance mutations. Check its inclusion rules, since novelty relative to earlier datasets is not the same as first report.
- **Option B:** a cutoff year chosen only from first-report-year counts so that at least 30 mutations fall after it.
- A fixed 2020 date is not allowed unless the counts support it.

### 1.6 Claims matrix with kill criteria

| Claim | Falsifiable hypothesis | Data condition | Kill criterion (proposed) | Fallback |
| :--- | :--- | :--- | :--- | :--- |
| C1 Benchmark | No HSV resistance benchmark with polymorphism negatives and grouped splits exists | Tiered labels | Prior art found | Reframe as an extension |
| C2 Leakage | Random splits inflate HSV performance vs position- and clade-grouped splits | Positions plus clade clusters | Gap under 0.03 AUROC with CI including 0 | Report as null |
| C3 Predictor | Beats zero-shot, GEMME and position-prior baselines at held-out positions | At least 15 resistant positions per gene-drug pair | Fewer than 15, or paired CI includes 0 | Pooled only |
| C4 Pathway gate | Routing helps most at low label fractions | At least 200 labeled rows; at least 3 independent sources | Interval includes 0 at 10 and 25 percent labels | Negative result |
| C5 Temporal | Pre-cutoff models rank later mutations | At least 30 post-cutoff mutations | Fewer than 30 | Descriptive |
| C6 Masking | TK loss-of-function rule beats additive and max-of-singles | At least 10 engineered doubles | Fewer than 10 | Case series |
| C7 Assay | Top-ranked UL5/UL52 predictions beat matched controls | Collaborator | No collaborator | Predictions only |
| C8 Derived outputs | n/a | Observed and predicted cells separate | Untested shown as negative | n/a |
| C9 Exploratory doubles and structure | n/a | Same search as C6 | n/a | Drop |
| C10 Multi-site head | Deferred | More than 30 doubles, at least 3 sources | n/a | Paper 3 |

### 1.7 Conflict rule (HerpesDRG, CHARMD, CMI Table S3, paper)
1. Prefer the higher evidence tier. 2. Prefer the original paper over a database summary. 3. If still tied, mark **CONFLICTED**, exclude from primary training and evaluation, and rerun with the row treated both ways. 4. Record both sources and the reason; report the number of conflicts. 5. Never default to "resistant."

### 1.8 Other frozen rules
- Untested mutation-drug pairs are absent, never negative.
- Co-occurrence in an isolate is never an "engineered double."
- Rows with unresolved curator notes are quarantined until a human reviews them.
- Doubles search: Day 8 is a checkpoint. If under 10 are found and the logged query list is exhausted, report "none found in the logged search," not "none exist."
- No "to our knowledge" statement without a logged search.
- Seal the held-out test positions by assigning them with a seeded script and hashing the assignment before modeling.

### 1.9 Standing Manuscript Integrity Rule (The Log-Traceability Invariant)
**Nothing enters the manuscript unless a log file contains it.**
- Every single number, point estimate, confidence interval, standard deviation, sample size, or metric reported in the manuscript text, abstract, tables, or supplementary materials **MUST trace directly to an existing, verifiable log file on disk** (`logs/*.json`, `logs/*.tsv`, `logs/*.log`).
- No metric or interval may be computed on-the-fly in prose, back-derived, or asserted from expectation.
- Every citation establishing an assay or IC50 cutoff must be verified against the primary published text before inclusion. Fabricated or proxy citations (such as citing a genotypic review for unverified phenotypic IC50 thresholds) are strictly prohibited.

---

## 2. Stage 1: Rights and coverage

| Check | Finding so far | Action |
| :--- | :--- | :--- |
| HerpesDRG license | MIT in the repository | Done; record commit hash |
| HerpesDRG gene coverage | Active HSV-1 rows: UL23 377, UL30 341, UL5 26, UL52 9, UL29 1. HSV-2: UL30 99, UL23 93, UL5 5. **UL8 absent** | Done |
| HerpesDRG status field | The database README says to use only status "A" | Done |
| CMI 2025 review Table S3 | Review exists; supplement cited. A corrigendum exists | Read the corrigendum, supplement license and redistribution terms |
| CHARMD | Includes both resistance mutations and natural polymorphisms. HSV coverage and terms unverified | Check HSV coverage and terms of use |
| Sauerbrei database | Exact citation unverified | Confirm citation and access |

If any source's terms are restrictive, release only identifiers, derived labels and scripts.

---

## 3. Stage 2: Extraction and source register

| Source | Provides | Notes |
| :--- | :--- | :--- |
| HerpesDRG-DB (github.com/ojcharles/herpesdrg-db) | UL23, UL30 labels; sparse UL5/UL52 | Hash the raw file. Record `control_strain` for clade grouping |
| CMI 2025 review, Table S3 | 1,604 mutation-drug combinations, including helicase-primase mutations for HSV-1 and HSV-2, and clinical correlation for 513 mutations. | Most important for UL5/UL52/PTV/AMNV |
| CHARMD | Curated resistance and polymorphism calls | After the rights check |
| Sauerbrei database | Historical baseline | Verification source |
| BAC study of AMNV resistance | Recombinant single-mutation UL5 and UL52 viruses. | Tier A for AMNV |
| Pritelivir combination study | K356N and K356Q with greater than 30-fold PTV and AMNV resistance. | Tier A or B check |
| Recombinant HSV-1 TK mutation study | Evaluates TK mutations (and doubles) in recombinant viruses. | Doubles search |
| UL29 P597L paper | A recombinant HSV-1 carrying P597L ICP8 was generated and tested. | Confirm PMID |
| 302-isolate HSV-1 survey | Sequence analysis of TK and polymerase from over 300 clinical isolates. | Polymorphism source; check whether isolates were phenotyped |
| UW 2026 clinical testing preprint | Six years of clinical genotypic testing with novel mutations. | Preprint; use for temporal analysis with care |
| NCBI Virus, RefSeq NC_001806.2 and NC_001798.2, UniProt, PDB, AlphaFold, InterPro | Sequences, structures, domains | Sequence is context, never a label |

**Rules:** no mutation, fold-change, PMID or citation from model memory. Save raw files with SHA-256. A citation is valid only if its DOI or PMID resolves and the title matches. Log every failure.

**Required fields on every row:** gene, virus type, position, ref and alt residue, `control_strain` or wild-type background, drug, raw value, `fold_change`, `threshold_used` and its source, assay type, `evidence_tier`, `phenotype_label`, `co_mutations_present`, `first_report_year`, source ID, PMID or DOI, curator note, extractor initials.

**Doubles search** (literature only): PubMed query for recombinant or marker-transfer viruses carrying both a TK and a polymerase mutation, and within UL5/UL52, with a measured joint phenotype. Log every query.

---

## 4. Stage 3: Quality control
- Deduplicate on (gene, position, substitution, drug, source), and remove the same isolate reported by several papers.
- A second person independently tiers about 50 random rows; report Cohen's kappa (proposed target above 0.80).
- Spot-check at least 10 percent of rows against the original papers.
- Maintain the conflict log (Section 1.7).
- Review all quarantined rows.

---

## 5. Stage 4: Counts and power

### 5.1 Baseline from the real HerpesDRG extract (v0, unverified, no thresholds applied)

Commit 22fc18676a1ea49cd6d7cd1c3280db210a7580e6, SHA-256 e1d9f61caa37f47b6a102ff020a4f483b2a82dc5735bd6688cda57827c36f770. 2,035 rows in total, 1,696 active, 949 active HSV mutations, 2,110 tested HSV mutation-drug rows.

| Gene, drug | Rows | Distinct positions | Positions labeled resistant | Positions labeled polymorphism | Fold-change rows |
| :--- | :--- | :--- | :--- | :--- | :--- |
| UL23, ACV | 468 | 212 | 76 | 105 | 163 |
| UL30, ACV | 432 | 305 | 17 | 208 | 161 |
| UL30, FOS | 402 | 291 | 12 | 218 | 133 |
| UL5, PTV | 20 | 6 | 0 | 1 | 19 |
| UL5, AMNV | 14 | 8 | 0 | 0 | 14 |
| UL52, PTV | 3 | 1 | 0 | 0 | 3 |
| UL52, AMNV | 6 | 5 | 0 | 0 | 6 |
| UL29, ACV | 1 | 1 | 0 | 0 | 1 |

"Labeled resistant" is the database's own call; thresholds applied later will raise these. Active HSV-1 UL23 unique mutations: 255 (292 across all statuses). No engineered doubles identified: the 8 co-mutation rows are isolates or lab-induced strains. First-report year, tier, and C1/C2 are not in the file.

### 5.2 Tables to produce
- **A:** per gene-drug pair: resistant positions (A+B), resistant positions (A only), polymorphism positions (C1, C2), fold-change rows, independent sources.
- **B:** label purity by position; positions with at least two substitutions of differing labels.
- **C:** engineered doubles by type; co-occurrence counted separately.
- **D:** labs, assay types and threshold variation per drug.
- **E:** mutations by first-report year and the count after the cutoff.
- **F:** HSV-1 versus HSV-2.

### 5.3 Strict additions
- **Threshold sensitivity:** recount resistant positions at several cited thresholds.
- **Source concentration:** top-source share and leave-one-source-out counts (one source supplies roughly 18 to 19 percent of UL30 rows).
- **Power check:** with about 15 positive positions at AUROC 0.8, a 95 percent interval is roughly ±0.14 (Hanley-McNeil), so use paired tests on grouped folds and claim no tiny differences.
- **Negative controls:** TK should not affect FOS, PTV or AMNV, and UL5/UL52 should not affect ACV. The extract has **2 positions labeled resistant for UL23 and FOS**; investigate each as a possible annotation error or missed co-mutation.

---

## 6. Stage 5: Feasibility probes (development folds only)
Fit trivial baselines: position prior, BLOSUM62, Grantham and conservation. If they approach model performance, the task is partly position memorization and C3 is reframed. The remaining gap is the realistic headroom.

---

## 7. Stage 6: Decision rules (frozen)

| Condition | Decision |
| :--- | :--- |
| Under about 15 resistant positions for a gene-drug pair | Pooled only; no per-drug claim |
| Under about 200 labeled rows overall | Drop C4 |
| Under about 30 post-cutoff mutations | C5 descriptive |
| High fraction of single-label positions | Position-shortcut baseline and within-position ranking mandatory |
| Fewer than 10 engineered doubles | C6 case series, no AUC |
| 10 to 30 doubles | C6 exploratory with exact intervals |
| More than 30 doubles from at least 3 sources | Consider the multi-site head |
| UL5/UL52 below minimum | Helicase-primase secondary |
| C2 rows exceed half the negatives | Sensitivity analysis excluding C2 |
| Conflicts exceed about 10 percent | Report prominently and temper claims |

---

## 8. Stage 7: Novelty audit

**Rule:** try to disprove each claim. Log date, database (PubMed, bioRxiv, arXiv, Semantic Scholar, Google Scholar, plus GitHub, Hugging Face, Zenodo and Kaggle), exact query, hits screened, the ten closest results with links, closest prior work, and a decision (**novel for HSV, novel method, not novel**). A second person repeats a sample of queries, a herpesvirology expert reviews, and the search is re-run before submission. Include HCMV and VZV phrasings and forward citations of PRIME, HerpesDRG and CHARMD.

### 8.1 Findings so far (first pass, not systematic; absence in my results is not proof)

| Claim | Closest prior work found | Preliminary read |
| :--- | :--- | :--- |
| C1 Benchmark | HerpesDRG is a database and annotation tool; CHARMD curates resistance and polymorphisms. No HSV benchmark with grouped splits seen | Likely novel for HSV; unconfirmed |
| C2 Leakage | PRIME (SARS-CoV-2): random-split R² above 0.90, about 0.23 on unseen sites, negative for base models. Regression on a large dataset | Concept exists; HSV-specific, clade-aware, multi-regime version appears open |
| C3 Predictor | No HSV predictor found. AMRscope predicts resistance from PLM embeddings of single mutations in bacteria. An MIT thesis showed language models can predict resistance in HIV-1 and beta-lactamase. | Method not novel; HSV application plausible |
| C4 Pathway gate | A drug-conditioned resistance model (LOGICA) reports modest held-out-gene gains. Gene-drug masking not yet searched | Pending; decides the claim |
| Contrastive plus LoRA | LoRA with contrastive objectives on viral proteins and supervised contrastive learning for resistance mutations exist | Not a method novelty |
| C5 Temporal | Not yet searched | Pending |
| C6 Masking | Not yet searched | Pending |
| Data sources | No deep mutational scanning for HSV targets found | Confirms sparsity |
| Boundary | HCMV, VZV resistance prediction not searched | Pending |

### 8.2 Searches still to run
1. Temporal or time-split evaluation in resistance mutation ranking (HIV, TB, influenza, SARS-CoV-2).
2. Gene-drug masked or routed models in multi-drug resistance prediction.
3. Masking epistasis or double-mutant models in resistance prediction.
4. HSV, HCMV, VZV resistance machine learning, protein language model, and contrastive terms.
5. GitHub, Hugging Face, Zenodo and Kaggle for HSV resistance datasets or models.
6. Forward citations of the HerpesDRG, PRIME and CHARMD papers.

---

## 9. Stage 8: Fact verification (source and PMID needed)

| Fact | Status |
| :--- | :--- |
| FCV is famciclovir; foscarnet is FOS | Confirm in a pharmacology reference |
| Pathway table (ACV needs TK; FOS, PTV, AMNV do not) | Cite sources. Penciclovir-resistant HSV-1 clones lacked polymerase mutations, consistent with TK dependence |
| UL8 role | Check whether any UL8 mutation was phenotyped |
| UL29 P597L | Confirm PMID and drug specificity |
| Polymorphism examples (UL30 E226D, T566A, A351D) | Verify each from a primary source |
| HerpesDRG counts | Use the current file, not the paper's older figures |
| ProtBERT context length (about 512 tokens) | Confirm; UL5, UL30 and UL52 need windowing |
| Corruption rate | PRIME used 15 percent masking; 2.5 percent is unsourced. Drop it or justify by your own ablation |
| LoRA trainable fraction | Measure it; do not assert a percentage |
| Context sequences | HerpesDRG maps mutations to the NCBI reference, so score against the reference by default and use isolate sequences only where available |
| Grant scope | Confirm with the PI that antiviral-resistance prediction fits the approved objectives |

---

## 10. Stage 9: Outputs

`AUDIT_REPORT.md` with: counts Tables A to F with raw snapshot hashes; the decision matrix (each claim Primary, Exploratory, Deferred or Dropped); the novelty log with decisions; the verified-facts table; kappa and the conflict log; pilot baseline results; and a one-paragraph locked scope statement.

**Locked scope template:**
> This paper presents an evidence-tiered, position- and clade-stratified benchmark for HSV resistance mutations (C1, C2). A drug-conditioned predictor is evaluated for [gene-drug pairs that pass the position threshold] (C3). The pathway-constrained gate is [Primary / Exploratory / Dropped] (C4). The temporal holdout is [Primary / Descriptive] (C5). The masking analysis is [Exploratory / Case series / Dropped] (C6). Cross-resistance and landscape outputs are derived predictions (C8). Multi-site epistasis modeling is deferred (C10). Every claim traces to the counts and the novelty log.

---

## 11. Gates, timeline and division of labor

| Gate | Condition |
| :--- | :--- |
| 1 Rights and coverage | Sources for UL5/UL52 identified, or helicase-primase demoted |
| 2 Data | Per-pair claims only where thresholds are met |
| 3 Feasibility | Baselines far below the model, or C3 reframed |
| 4 Novelty | Every abstract claim has a logged decision |
| 5 Sign-off | PI and second annotator approve the locked scope |

| Days | Work |
| :--- | :--- |
| 1 | Freeze Stage 0; rights and coverage |
| 2 to 8 | Extraction and QC (doubles checkpoint on Day 8) |
| 3 to 10 | Novelty searches in parallel |
| 9 to 11 | Counts, power, feasibility probes, decision rules |
| 12 to 14 | Fact verification and the report |
| 15 and later | Modeling only after sign-off |

**Agent (mechanical):** extraction scripts, downloads with hashes, DOI and first-year lookups, deduplication, counts and tests, running and logging queries.  
**Human:** tiering, thresholds, C1 versus C2, conflict resolution, quarantined rows, novelty decisions, fact verification, spot-checks (at least 10 percent), the locked scope.

---

## 12. Known open issues
- The thresholds, temporal rule choice, and HSV-2 scope still have to be decided by you and the PI.
- UL5/UL52 and engineered-double data are likely thin; the audit may demote helicase-primase and C6.
- CHARMD coverage and the CMI supplement's license and corrigendum are unchecked.
- Every novelty search is still preliminary.
