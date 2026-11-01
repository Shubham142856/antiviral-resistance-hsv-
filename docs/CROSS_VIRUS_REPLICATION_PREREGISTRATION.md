# Pre-Registration: Cross-Herpesvirus Replication of the Target-Specific Radicalness Dissociation

**Document Version:** 1.0 (Frozen Pre-Registration)  
**Date:** 2026-11-02  
**Target:** Cross-Herpesvirus Replication (HCMV UL97/UL54 and VZV ORF36/ORF28)  
**Status:** Pre-registered BEFORE running the replication AUROC evaluation.  
**Repository:** [`https://github.com/Shubham142856/antiviral-resistance-hsv-`](https://github.com/Shubham142856/antiviral-resistance-hsv-)

---

## 1. Rationale & Pre-Registered Predictions

In HSV-1, amino acid radicalness (BLOSUM62) predicts resistance on the dispensable prodrug activator (thymidine kinase UL23, AUROC ~ 0.6130–0.6380) but operates at chance on the essential replicative enzyme (DNA polymerase UL30, AUROC ~ 0.4531–0.5074).

We evaluate whether this target-dependent rank-AUROC dissociation replicates across orthologous systems in human cytomegalovirus (HCMV) and varicella-zoster virus (VZV):
1. **HCMV (Human Betaherpesvirus 5):**
   - **UL97 (Phosphotransferase / Prodrug Activator):** Functionally analogous to HSV-1 UL23 in prodrug phosphorylation. Resistant variants frequently arise via catalytic loss-of-function.
     - **Pre-Registered Prediction:** BLOSUM62 radicalness **clears chance** ($\text{AUROC} > 0.55$, direction: higher radicalness -> higher resistance probability).
   - **UL54 (DNA Polymerase Catalytic Subunit):** Essential replicative polymerase, homologous to HSV-1 UL30.
     - **Pre-Registered Prediction:** BLOSUM62 radicalness **operates near chance** ($\text{AUROC} \in [0.45, 0.55]$).
2. **VZV (Human Alphaherpesvirus 3):**
   - **ORF36 (Thymidine Kinase):** Functionally and structurally homologous to HSV-1 UL23 ($N_{\text{tot}}=28, N_{\text{res}}=18$).
     - **Pre-Registered Prediction:** BLOSUM62 radicalness **clears chance** ($\text{AUROC} > 0.55$).
   - **ORF28 (DNA Polymerase):** Homologous to HSV-1 UL30. **Power Caveat:** With only $N_{\text{res}}=4$ resistant missense variants in HerpesDRG-DB, ORF28 is underpowered; its AUROC is documented as descriptive context only.

---

## 2. Dataset Extraction & Inclusion Criteria

- **Data Source:** `herpesdrg-db/herpesdrg-db.tsv` (authoritative curated database).
- **Inclusion Rules:**
  - Strictly missense substitutions with parseable wildtype, position, and mutant residues.
  - Deduplicated on unique amino acid substitutions `(gene, aa_change)`.
  - Resistance threshold for HCMV: fold-change $\ge 2.0\times$ (standard clinical cut-off for ganciclovir/foscarnet/cidofovir per Limoges/ANRS and HerpesDRG) or reported 'Resistant'. Sensitive/polymorphism if fold-change $< 2.0\times$ or reported 'Polymorphism'/'Sensitive'.
  - Resistance threshold for VZV: fold-change $\ge 3.0\times$ or reported 'Resistant'.

---

## 3. Explicit Falsification Criteria

- **Kinase Dissociation Falsifier:** If BLOSUM62 radicalness on HCMV UL97-GCV yields $\text{AUROC} \le 0.50$, the hypothesis that prodrug activators tolerate and reflect non-specific radical disruptions is **falsified**.
- **Polymerase Dissociation Falsifier:** If BLOSUM62 radicalness on HCMV UL54 (GCV, FOS, or CDV) yields $\text{AUROC} \ge 0.60$ with 95% bootstrap CI strictly excluding 0.50, the claim that essential replicative polymerases selectively restrict radical disruptions to chance is **falsified**.

---

## 4. Execution Commitment

The evaluation script `scripts/replication/evaluate_cross_virus_replication.py` will be executed immediately following the recording of this protocol.
