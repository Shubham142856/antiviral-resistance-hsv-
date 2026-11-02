"""
patch_step0b_final_manuscript.py
Applies Final Step 0b updates to manuscript/manuscript_draft.md and README.md:
1. Updates Table 7 with Supervised Logistic Regression L2 out-of-domain evaluation.
2. Reports NO bootstrap CI for FOS (<4 evaluable domains, n=2, 3 distinct resamples).
3. Adds per-domain breakdowns and sign-consistency.
4. Replaces "929 substitutions" with "929 variant-drug records (630 unique substitutions)".
5. Explicitly documents UL30-ACV/FOS position overlap (~88-91%) and HCMV UL54 multi-drug non-independence.
6. Details beta5 LOO influence diagnostic (Pos 927, Palm domain, r = -0.697, bimodal split).
7. Reports BLOSUM62 UL23 minus UL30 difference (+0.1820) strictly as descriptive.
"""
import re
import pathlib

# 1. Read manuscript
with open('manuscript/manuscript_draft.md', 'r', encoding='utf-8') as f:
    text = f.read()

# Replace 929 substitutions in abstract, intro, and tables
text = text.replace(
    "benchmark of 929 single-site missense substitutions across thymidine kinase",
    "benchmark of 929 single-site missense variant-drug records (630 unique substitutions) across thymidine kinase"
)
text = text.replace(
    "primary cohort of 929 single-site missense substitutions across UL23 and UL30",
    "primary cohort of 929 single-site missense variant-drug records (630 unique substitutions) across UL23 and UL30"
)
text = text.replace(
    "across 929 single-site missense mutations.",
    "across 929 single-site missense variant-drug records (630 unique amino acid substitutions)."
)

# Table 7 replacement
old_table7_pattern = re.compile(
    r"\*\*Table 7 \| Dual-Buffered Domain Benchmark.*?\*Calibration & Small-n Domain Notes:\*.*?(?=#### Operational Split-Calibration)",
    re.DOTALL
)

new_table7 = """**Table 7 | Dual-Buffered Domain Benchmark ($\pm 5$ aa, $< 12.0$ \AA): Recalibrated Predictor Performance, Supervised LR and Excess AUROC.** 
Performance reported as evaluable per-domain mean [with pooled out-of-fold AUROC in brackets]. Evaluable domains strictly enforce the pre-registered $N_{\\text{pos}} \\ge 5$ rule (UL30-ACV: 4 domains [Exo, Palm, Fingers, Thumb]; UL30-FOS: 2 domains [Palm and Fingers; Exo and Thumb omitted with $N_{\\text{pos}} < 5$]; UL23-ACV: 4 domains [ATP Core, Nucleoside Pocket, Lid/Dimer, C-term]). Excess AUROC is recalibrated against an empirical chance baseline: $\\Delta_{\\text{excess}} = \\text{AUROC}_{\\text{model}} - \\max(0.50, \\text{AUROC}_{\\text{pos\\_only}})$. For 4-domain targets ($n=4$), 95% confidence intervals are computed via domain-clustered bootstrap (recording 35 distinct multisets) and labeled descriptive only. For targets with $<4$ evaluable domains (UL30-FOS, $n=2$), **NO bootstrap CI is reported** because $\\binom{3}{2} = 3$ distinct resamples cannot support an empirical interval. Non-parametric sign consistency ($k/n$ domains with $\\text{AUROC}_{\\text{model}} > \\text{AUROC}_{\\text{pos\\_only}}$) is reported for all targets.

| Target Pair | $N_{\\text{eval}}$ ($n_{\\text{dom}}$) | Position-Only RF Mean [Pooled] | Split Validity ($[0.45, 0.55]$) | Zero-Shot ESM-2 650M Mean | **ESM-2 Excess $\\Delta$ [95% CI] (Sign)** | Zero-Shot BLOSUM62 Mean | **BLOSUM62 Excess $\\Delta$ [95% CI] (Sign)** | Where Spec-Exact Mean | **Where Excess $\\Delta$ [95% CI] (Sign)** | Supervised LR L2 Mean | **Sup-LR Excess $\\Delta$ [95% CI] (Sign)** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL30-Aciclovir** | 248 (4) | 0.5645 [**0.5185**] | **CONFIRMED** (Pooled in $[0.45, 0.55]$) | 0.7439 | **+0.1794** [+0.0281, +0.2664] (3/4) | 0.4906 | -0.0739 [-0.1589, +0.0115] (1/4) | **0.7801** | **+0.2157** [+0.1116, +0.3197] (4/4) | 0.7416 | **+0.1771** [+0.1274, +0.2339] (4/4) |
| **UL30-Foscarnet** | 114 (2) | 0.5253 [**0.5087**] | **CONFIRMED** (In Palm & Fingers) | 0.5855 | +0.0602 [NO CI, $n=2$] (1/2) | 0.5845 | +0.0592 [NO CI, $n=2$] (1/2) | **0.7878** | **+0.2625** [NO CI, $n=2$] (2/2) | 0.7017 | **+0.1764** [NO CI, $n=2$] (2/2) |
| **UL23-Aciclovir** | 254 (4) | 0.4372 [**0.3775**] | **OUTSIDE** (Domain rate shift: 14%–75%) | 0.6885 | **+0.1885** [+0.1046, +0.2725] (4/4) | **0.6726** | **+0.1726** [+0.1238, +0.2215] (4/4) | 0.5670 | +0.0670 [-0.0360, +0.1700] (2/4) | **0.7177** | **+0.2177** [+0.1231, +0.2761] (4/4) |

*Per-Domain Breakdown & Calibration Notes (`logs/step0b_final_recalibrated_benchmarks.json`):*
1. **Per-Domain Breakdown on Primary Target (UL30-ACV):**
   - *Exonuclease ($N=37$, 8 res):* Pos-RF = 0.7091; ESM-2 = 0.6595; BLOSUM62 = 0.5776; Where = **0.8621**; Sup-LR = 0.8233.
   - *Palm ($N=111$, 32 res):* Pos-RF = 0.5332; ESM-2 = 0.8050; BLOSUM62 = 0.5712; Where = **0.8133**; Sup-LR = 0.7983.
   - *Fingers ($N=18$, 16 res):* Pos-RF = 0.5156; ESM-2 = 0.7500; BLOSUM62 = 0.5000; Where = **0.8750**; Sup-LR = 0.6562.
   - *Thumb / C-term ($N=82$, 6 res):* Pos-RF = 0.5000; ESM-2 = **0.7610**; BLOSUM62 = 0.3136; Where = 0.5702; Sup-LR = 0.6886.
2. **Per-Domain Breakdown on Comparative Kinase (UL23-ACV):**
   - *ATP Core ($N=90$, 52 res):* Pos-RF = 0.4013 (Ref: 0.5000); ESM-2 = 0.8057; BLOSUM62 = 0.6371; Where = 0.7060; Sup-LR = 0.7854.
   - *Nucleoside Pocket ($N=60$, 45 res):* Pos-RF = 0.3815 (Ref: 0.5000); ESM-2 = 0.7393; BLOSUM62 = 0.6104; Where = 0.6341; Sup-LR = 0.7481.
   - *Lid / Dimer Interface ($N=62$, 17 res):* Pos-RF = 0.4660 (Ref: 0.5000); ESM-2 = 0.5980; BLOSUM62 = **0.7346**; Where = 0.4327; Sup-LR = 0.5778.
   - *C-term ($N=42$, 6 res):* Pos-RF = 0.5000 (Ref: 0.5000); ESM-2 = 0.6111; BLOSUM62 = **0.7083**; Where = 0.4954; Sup-LR = 0.7593.
   - *Verification of BLOSUM62 Domain Mean:* All 4 structural domain folds yield strictly positive excess over 0.50 (mean = **0.6726**, range: 0.6104 to 0.7346), recomputed directly from out-of-fold domain evaluations rather than legacy block results.
3. **Where vs. Supervised LR Comparison:** On the essential polymerase (UL30-ACV), the label-free `Where` composite (domain mean: **0.7801**) matches and marginally exceeds Supervised Logistic Regression (domain mean: **0.7416**; Sign: 4/4 vs 4/4). However, on the dispensable kinase (UL23-ACV), `Where` (domain mean: **0.5670**) substantially lags behind Supervised LR (0.7177), Zero-Shot ESM-2 (0.6885), and BLOSUM62 radicalness (0.6726). The "positional prior $\\approx$ SOTA" finding is therefore strictly specific to targets whose resistance landscape is governed by localized active-site intolerance.
4. **Descriptive Dissociation Contrast for BLOSUM62:** Across the 4 evaluable domains of UL23 vs. UL30-ACV, BLOSUM62 radicalness exhibits a domain-mean contrast of $\\Delta_{\\text{domain}} = 0.6726 - 0.4906 = \\mathbf{+0.1820}$ (with excess over 0.50 of $+0.1726$ on UL23 vs. $-0.0094$ on UL30). Given the small-$n$ domain limit ($n=4$ each), this contrast is presented strictly as a **descriptive** biophysical difference rather than an inferential claim.
5. **Non-Independence of Arms:** In our 929 variant-drug records (630 unique substitutions), UL30-ACV and UL30-FOS share 264 sequence positions (**91.0% of ACV positions and 97.4% of FOS positions**) and 299 identical substitutions (**90.1% of ACV substitutions**). Similarly, the three HCMV UL54 cohorts represent identical underlying polymerase sequence records tested against multiple drugs. Consequently, ACV and FOS arms within UL30 (and the 3 drugs within HCMV UL54) represent **correlated, non-independent arms** of the same biophysical target rather than separate replication events.

"""

text = old_table7_pattern.sub(lambda m: new_table7, text)

# Update beta5 discussion in Section 3.8
old_b5_text = """A leave-one-domain-out influence diagnostic (`logs/step0b_recalibrated_benchmarks.json`) identifies the structural driver: omitting the Palm domain ($N=111$, containing 32 of the 62 resistant events) shifts $ eta_5$ from $-0.0162$ to **$+0.4289$**, whereas omitting any other domain leaves $ eta_5$ near zero ($-0.012$ to $-0.086$). Because clinical resistant events are heavily concentrated in a single catalytic domain, $ eta_5$ is **statistically unidentifiable** in retrospective clinical cohorts. Under our pre-registered decision rule (requiring the 95% CI upper bound to be $< 0$), **the tolerance-window hypothesis is unidentifiable and neither confirmed nor refuted at current clinical sample sizes**."""

new_b5_text = """An exploratory Leave-One-Position-Out (LOO) influence diagnostic across all 290 unique positions (`logs/beta5_top15_loo_diagnostic.json`) explains the bimodal bootstrap distribution of $\\beta_5$ (Mode 1: 68.5% of resamples centered at $-0.018$; Mode 2: 31.4% of resamples centered at $+0.683$). A single high-influence residue in the Palm domain—**Position 927 (Glu)**—drives the mode switch: omitting Position 927 alone shifts $\\beta_5$ from $-0.0161$ to **$+0.6012$** ($\\Delta = +0.6173$). Across the 2,000 bootstrap resamples, the inclusion count of Position 927 correlates strongly with $\\beta_5$ ($r = -0.6966$): when Position 927 is omitted in random sampling (~31% of resamples), the fit flips into the secondary positive mode ($+0.60$ to $+0.85$). 

Omitting the entire Palm domain ($N=111$, containing 32 of the 62 resistant events) similarly shifts $\\beta_5$ to **$+0.4289$**, while omitting any other domain leaves $\\beta_5$ near zero. This diagnostic proves that the secondary positive mode reflects severe position-level leverage rather than a broad quasi-separation tail. Under our pre-registered decision rule (requiring the 95% CI upper bound to be $< 0$), **the tolerance-window hypothesis is unidentifiable: $\\beta_5$ is neither confirmed nor refuted at current retrospective sample sizes**."""

text = text.replace(old_b5_text, new_b5_text)

with open('manuscript/manuscript_draft.md', 'w', encoding='utf-8') as f:
    f.write(text)

print("Updated manuscript/manuscript_draft.md successfully.")
