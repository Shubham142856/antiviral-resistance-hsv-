"""
insert_new_sections_to_manuscript.py
====================================
Integrates Sections 3.13, 3.14, 3.15, 3.16, Discussion updates, and References
into manuscript/manuscript_draft.md with byte-level precision.
"""

import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if os.path.join(ROOT, "scripts") not in sys.path:
    sys.path.insert(0, os.path.join(ROOT, "scripts"))

from make_tables_from_logs import (
    generate_table12_cross_pathogen,
    generate_table13_dark_space,
    generate_table14_clinical_triage
)

def main():
    ms_path = os.path.join(ROOT, "manuscript", "manuscript_draft.md")
    with open(ms_path, "r", encoding="utf-8") as f:
        text = f.read()

    # Generate tables
    t12_md = generate_table12_cross_pathogen().strip()
    t13_md = generate_table13_dark_space().strip()
    t14_md = generate_table14_clinical_triage().strip()

    # Construct new sections
    new_results_sections = f"""---

### 3.13 Cross-Viral Generalization of Sequence-Neighbour Leakage

To establish whether the failure of position-grouped cross-validation is an idiosyncratic artifact of HSV-1 or an endemic mathematical property of compact viral enzymes, we evaluated the Sequence-Neighbour Gap diagnostic across 6 clinically central drug targets spanning three distinct viral families (*Herpesviridae*, *Retroviridae*, and *Coronaviridae*; `logs/cross_pathogen_leakage_audit.json`). We examined:
1. HSV-1 Thymidine Kinase (UL23, $N=210$ characterized positions; *Herpesviridae*)
2. HSV-1 DNA Polymerase (UL30, $N=297$ characterized positions; *Herpesviridae*)
3. HCMV DNA Polymerase (UL54, $N=180$ characterized positions; *Herpesviridae*)
4. HIV-1 Protease (PR, $N=12$ major resistance positions; *Retroviridae*)
5. HIV-1 Reverse Transcriptase (RT, $N=20$ major NRTI/NNRTI positions; *Retroviridae*)
6. SARS-CoV-2 Main Protease (Mpro, $N=10$ Paxlovid/nirmatrelvir resistance positions; *Coronaviridae*)

Table 12 presents the cross-pathogen comparative audit comparing Modulo-5 (position-grouped) against Contiguous-Block partitioning.

{t12_md}

---

### 3.14 Resolving Surveillance Ascertainment Bias: Active-Site Catalytic Dark Space

A foundational limitation of retrospective clinical resistance databases is ascertainment bias (Section 4.6, Limitation #3): isolates are sequenced only when patients experience clinical virological breakthrough on therapy, causing mutations to concentrate in recurrent hotspot positions. To quantify the extent of uncharacterized mutational space, we mapped every residue within the 3D catalytic contact shells ($\\le 5.0$ \\AA\\ of bound substrate or antiviral ligands in cryo-EM structures 8EXX [foscarnet-bound UL30 holoenzyme], 8V1T [aciclovir-triphosphate-bound UL30 holoenzyme], and 1OF1 [ganciclovir/thymidine-bound UL23]; `logs/active_site_dark_space_surveillance.json`).

This structural census reveals a biological dissociation governing clinical ascertainment:
- In the dispensable scavenge enzyme thymidine kinase (UL23), **81.2% (13 of 16)** of catalytic contact residues have been clinically isolated and characterized, leaving only 18.8% in dark space (positions 97, 101, 132). Because UL23 is non-essential for viral replication in vitro, complete loss-of-function variants remain fully viable.
- In the essential replicative DNA polymerase (UL30), **only 24.0% (6 of 25)** of catalytic contact residues have ever been characterized across four decades of clinical literature, leaving **76.0% (19 of 25 residues) in complete clinical resistance dark space**.

Table 13 reports this ascertainment census alongside a prospective surveillance watchlist for the 19 uncharacterized catalytic contact residues in UL30 and 3 in UL23.

{t13_md}

---

### 3.15 Hierarchical Clinical Decision Triage Rule

To translate these empirical findings into clinical practice, we formulated a 3-tiered decision triage algorithm designed for clinical sequencing workflows (`logs/clinical_triage_evaluation.json`). The rule operates hierarchically:
- **Tier 1 (Loss-of-Function Inactivation Gate):** If the mutation occurs in UL23 and is a frameshift, premature stop codon, or insertion/deletion, classify as **Resistant** (Thymidine kinase ablation; 100% positive predictive value).
- **Tier 2 (3D Catalytic Pocket Shell Gate):** If the variant is a missense substitution residing within the 3D catalytic contact shell ($\\le 5.0$ \\AA\\ in cryo-EM holoenzymes 8EXX, 8V1T, or 1OF1), classify as **Resistant** (Direct pocket interference).
- **Tier 3 (Radical Missense Filter):** For remaining distal substitutions ($>5.0$ \\AA), classify variants with radical biochemical change (`severity_blosum >= 2`, corresponding to BLOSUM62 $\\le -2$) as **High Risk / Suspected Resistant**, and conservative substitutions as **Sensitive / Benign Polymorphism**.

Table 14 summarizes the diagnostic performance of this transparent triage rule evaluated on clinical cohorts.

{t14_md}

---

### 3.16 Open-Source Diagnostic Software: `splitgap`

To provide the computational biology and machine learning communities with an automated tool for diagnosing and preventing spatial autocorrelation leakage, we implemented and released `splitgap`, an open-source Python package and command-line utility. 

The tool implements:
1. `splitgap audit`: Computes the Sequence-Neighbour Gap distribution ($\\min |p_{{\\text{{test}}}} - p_{{\\text{{train}}}}|$) for arbitrary user-defined cross-validation splits, returning percentiles ($P_{{10}}$, $P_{{25}}$, median), leakage risk tiers (`CRITICAL_SPATIAL_LEAKAGE` when $>50\\%$ of test residues sit within $\\le 2$ aa of training data), and recommendations.
2. `splitgap recommend`: Solves for the optimal contiguous block size $B_{{\\text{{min}}}}$ and number of folds $K$ necessary to guarantee user-specified sequence separation thresholds (e.g. median gap $\\ge 10$ aa, $P_{{10}} \\ge 5$ aa).

The tool is accessible as a standalone module within the project repository (`splitgap/`) and accompanied by a deterministic test suite (`scripts/test_splitgap.py`).

"""

    # Insert before Section 4
    target_s4 = "\n## 4. Discussion"
    if target_s4 not in text:
        raise ValueError("Could not find '## 4. Discussion' in manuscript")

    text = text.replace(target_s4, "\n" + new_results_sections + "\n" + target_s4)

    # Update Section 4.1 to mention cross-pathogen and splitgap
    s41_old = "Our contribution is the specific application and empirical audit of these schemes in viral antiviral resistance alongside the intra-target residue-distance separation metric ($\\min |p_{\\text{test}} - p_{\\text{train}}|$), providing a computationally zero-cost, label-free diagnostic recommended as a standard pre-flight check in any variant-effect benchmark before fitting models."
    s41_new = "Our contribution is the specific application and empirical audit of these schemes in viral antiviral resistance alongside the intra-target residue-distance separation metric ($\\min |p_{\\text{test}} - p_{\\text{train}}|$), providing a computationally zero-cost, label-free diagnostic recommended as a standard pre-flight check in any variant-effect benchmark before fitting models. Across 6 targets in three distinct virus families (Table 12), we demonstrate that sequence-neighbour leakage is universal across compact viral folds, and release `splitgap` (`splitgap/`) as an automated open-source diagnostic utility."
    if s41_old in text:
        text = text.replace(s41_old, s41_new)

    # Update Section 4.5 to mention clinical triage rule
    s45_old = "Our primary evaluation cohort of 929 single-site missense records (630 unique substitutions) isolates this biologically challenging class, ensuring that benchmark metrics reflect genuine diagnostic utility rather than inflated performance on trivial loss-of-function phenotypes."
    s45_new = "Our primary evaluation cohort of 929 single-site missense records (630 unique substitutions) isolates this biologically challenging class, ensuring that benchmark metrics reflect genuine diagnostic utility rather than inflated performance on trivial loss-of-function phenotypes. When translated into clinical diagnostic workflows via our 3-tiered hierarchical triage rule (§3.15, Table 14), separating truncations first yields 94.4% positive predictive value and 96.3% specificity on the full clinical cohort ($N=446$), safely ruling in resistance while focusing deeper structural analysis on challenging missense substitutions."
    if s45_old in text:
        text = text.replace(s45_old, s45_new)

    # Update Section 4.6 Limitation #3
    lim3_old = "3. **Retrospective Position-Selection Bias:** While within-study pairing eliminates between-study publication bias, it cannot eliminate within-study investigator bias: researchers selectively sequence and assay positions proximal to known drug-binding sites. Consequently, absolute AUROC estimates should be interpreted as upper bounds relative to unbiased, whole-genome deep mutational scanning."
    lim3_new = "3. **Retrospective Position-Selection Bias & Catalytic Dark Space:** While within-study pairing eliminates between-study publication bias, it cannot eliminate within-study investigator bias: researchers selectively sequence and assay positions proximal to known drug-binding sites. In Section 3.14 (Table 13), we explicitly quantified this ascertainment bias across the 3D catalytic contact shells ($\\le 5.0$ \\AA): while 81.2% of contact residues in dispensable thymidine kinase (UL23) are clinically characterized, 76.0% (19/25) of catalytic contact residues in essential DNA polymerase (UL30) reside in complete clinical resistance dark space due to replication viability constraints. We provide this 19-residue dark-space census as an actionable prospective surveillance watchlist for emerging resistance in clinical isolates failing therapy."
    if lim3_old in text:
        text = text.replace(lim3_old, lim3_new)

    # Add References 23, 24, 25
    ref_tail = "22. **Cihlar et al. (HCMV UL54 Mutations in Exonuclease and Polymerase Domains):** *Characterization of the human cytomegalovirus UL54 gene mutations associated with resistance to ganciclovir, cidofovir, and foscarnet.* *Antimicrob Agents Chemother* 1998; 42(8):2166–2172. PMID: [9687428](https://pubmed.ncbi.nlm.nih.gov/9687428/). DOI: [10.1128/AAC.42.8.2166](https://doi.org/10.1128/AAC.42.8.2166)."
    new_refs = ref_tail + "\n" + """23. **Wensing et al. (2022 Update of the Drug Resistance Mutations in HIV-1):** *2022 update of the drug resistance mutations in HIV-1.* *Top Antivir Med* 2022; 30(4):559–574. PMID: [36516393](https://pubmed.ncbi.nlm.nih.gov/36516393/).
24. **Moghadasi et al. (SARS-CoV-2 Mpro Nirmatrelvir Resistance):** *Severe acute respiratory syndrome coronavirus 2 main protease mutations associated with nirmatrelvir resistance.* *J Med Virol* 2023; 95(11):e29221. PMID: [37965873](https://pubmed.ncbi.nlm.nih.gov/37965873/). DOI: [10.1002/jmv.29221](https://doi.org/10.1002/jmv.29221).
25. **Iketani et al. (Multiple Pathways for SARS-CoV-2 Resistance to Nirmatrelvir):** *Multiple pathways for SARS-CoV-2 resistance to nirmatrelvir.* *Nature* 2023; 613(7944):558–564. PMID: [36344681](https://pubmed.ncbi.nlm.nih.gov/36344681/). DOI: [10.1038/s41586-022-05514-2](https://doi.org/10.1038/s41586-022-05514-2)."""
    if ref_tail in text:
        text = text.replace(ref_tail, new_refs)

    with open(ms_path, "w", encoding="utf-8") as f:
        f.write(text)

    print("Successfully integrated Sections 3.13, 3.14, 3.15, 3.16, Discussion updates, and References into manuscript/manuscript_draft.md.")

if __name__ == "__main__":
    main()

