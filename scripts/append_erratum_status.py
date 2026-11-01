"""
append_erratum_status.py
Appends Section 10 Erratum to STATUS.md
"""
import pathlib

p = pathlib.Path('STATUS.md')
erratum = """

---

## 10. Audit Erratum and Implementation Reconciliation Log (2026-11-02)

| Item | Focus Area | Action Taken | Evidence / Verification |
| :--- | :--- | :--- | :--- |
| **E10.1** | **UL23 Linear Disruption Sign Resolution (\\beta_2)** | In the initial committed results (`logs/where_what_test_results.json`), substitution severity was computed from raw masked marginals $\\log P(\\text{mut}) - \\log P(\\text{wt}) \\le 0$ without negation, making higher numerical values more benign and disruptive substitutions negative, resulting in $\\beta_2 = -0.3629$ [$-0.6831, -0.1235$]. Under the pre-registered convention where $\\Delta\\text{LLR} = \\log P(\\text{wt}) - \\log P(\\text{mut}) \\ge 0$ (higher = more disruptive), the correctly signed linear disruption coefficient on the dispensable kinase (UL23-ACV) is $\\mathbf{\\beta_2 = +0.3629}$ [$\\mathbf{+0.1235, +0.6831}$] ($p < 0.05$), passing the pre-registered directional prediction $\\beta_2 > 0$ with zero curvature ($\\beta_3 = +0.0423$ [$-0.1485, +0.2463$]). Both versions are recorded; no silent overwrites. | `logs/where_what_test_results.json`, `logs/where_what_spec_exact_results.json`, Git commit `273899b` |
| **E10.2** | **Where Feature Specification Deviation & Spec-Exact Rerun** | In the initial tolerance model run (`run_tolerance_models.py`), `where_composite` used raw $\\mu_{\\text{raw}}$ without negating it, effectively taking $-Z(\\text{Constraint}) - Z(d_{\\text{site}})$, which suppressed full model fit AUROC to 0.5993 and left $\\beta_1 = +0.1759$. When rerun with spec-exact Where ($Z(\\text{Constraint}) - Z(d_{\\text{site}})$ in `scripts/where_what/run_spec_exact_tolerance_models.py`):<br/>- Full model fit AUROC jumps to **0.8611** (UL30-ACV), **0.9168** (UL30-FOS), and **0.7551** (UL23-ACV).<br/>- $\\beta_1$ (Where) becomes massive and positive: $\\mathbf{+1.9372}$ [$\\mathbf{+1.0555, +2.7415}$] on UL30-ACV and $\\mathbf{+2.0624}$ [$\\mathbf{+1.5512, +2.7406}$] on UL30-FOS.<br/>- $\\beta_5$ on UL30-ACV is $\\mathbf{-0.0162}$ [$\\mathbf{-0.1876, +1.0402}$]. Since the 95% position-cluster CI spans zero, the tolerance-window hypothesis remains unconfirmed / not supported at this sample size ($N_{\\text{pos}}=62$). Both runs are archived and reported. | `logs/where_what_spec_exact_results.json` |
| **E10.3** | **Severity Axis Identifiability & Multiplicity Audit** | Audited position multiplicity across benchmark arms: On UL30-ACV, 86.9% of evaluated positions (252/290) possess exactly 1 variant, and only 2 positions (0.7%) have mixed resistant/sensitive labels (FOS: 235/271 single-variant, 3 mixed). On UL23, 60 positions (28.6%) have $\\ge 2$ variants (11 mixed). This establishes that the lack of within-position severity signal on UL30 is driven by clinical data scarcity (the severity axis is largely unidentifiable) rather than proving that substitution severity is biologically irrelevant. | `data/processed/primary_position_folds.tsv` |
| **E10.4** | **HCMV Radicalness Contrast Bootstrap Test** | Tested the empirical contrast between HCMV UL97 (GCV, Kinase) and UL54 (GCV, Polymerase) via 2,000 bootstrap resamples: Observed $\\Delta = \\mathbf{+0.1630}$, 95% bootstrap CI on $\\Delta$: $[\\mathbf{+0.0341, +0.2861}]$, empirical $p = 0.0070$ (strictly excludes zero). Replication claim calibrated to $n=1$ independent pair (HCMV); VZV polymerase arm ($N_{\\text{res}}=4$) is documented as not evaluable for contrast testing. | `herpesdrg-db/herpesdrg-db.tsv`, `logs/cross_virus_replication_results.json` |
"""

with open(p, 'a', encoding='utf-8') as f:
    f.write(erratum)
print("Section 10 Erratum successfully added to STATUS.md")

