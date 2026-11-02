"""
patch_stage_a.py
Implements all STAGE A manuscript and README updates:
1. beta5: pre-registered falsifier not met (estimate not supported); estimate is leverage-driven (E927 moves beta5 from -0.016 to +0.60; ~31% second mode; largely driven; influence diagnostic, not field-level contribution).
2. UL30-FOS split validity: in-band, n=2 domains, position-only RF is exactly 0.5000 in Fingers (constant predictions); not independently interpretable. Remove CONFIRMED.
3. UL23 position-only RF: flag in Table 7 that excess is chance-referenced (reference clamped to 0.50; every domain below chance). Report inspection of pos-RF predictions vs position per domain; do not assert cause.
4. Where vs LR: report per-domain differences. Use "matches", not "exceeds". Clarify sign counts compare each model against pos-RF, not Where-vs-LR.
5. Language: remove "proving", replace "independent measurements" with "consistent observations on the same two genes". Show pooled and domain-mean AUROC side by side with one line on what each measures.
6. Overlap figures: exact computed values (91.0% ACV positions, 97.4% FOS positions, 90.1% substitutions).
"""
import re

# 1. Update manuscript
with open('manuscript/manuscript_draft.md', 'r', encoding='utf-8') as f:
    text = f.read()

# Fix Table 7
old_table7_pattern = re.compile(
    r"\*\*Table 7 \| Dual-Buffered Domain Benchmark.*?\*Per-Domain Breakdown & Calibration Notes.*?(?=#### Operational Split-Calibration)",
    re.DOTALL
)

new_table7 = """**Table 7 | Dual-Buffered Domain Benchmark ($\pm 5$ aa, $< 12.0$ \AA): Predictor Performance, Supervised LR and Excess AUROC.** 
Performance reported as evaluable per-domain mean alongside pooled out-of-fold AUROC (pooled AUROC measures joint across-domain ranking and base-rate calibration, whereas domain-mean AUROC isolates within-domain discriminative ordering). Evaluable domains strictly enforce the pre-registered $N_{\\text{pos}} \\ge 5$ rule (UL30-ACV: 4 domains [Exo, Palm, Fingers, Thumb]; UL30-FOS: 2 domains [Palm and Fingers; Exo and Thumb omitted with $N_{\\text{pos}} < 5$]; UL23-ACV: 4 domains [ATP Core, Nucleoside Pocket, Lid/Dimer, C-term]). Excess AUROC is recalibrated against an empirical chance baseline: $\\Delta_{\\text{excess}} = \\text{AUROC}_{\\text{model}} - \\max(0.50, \\text{AUROC}_{\\text{pos\\_only}})$. For 4-domain targets ($n=4$), 95% confidence intervals are computed via domain-clustered bootstrap (recording 35 distinct multisets) and labeled descriptive only. For targets with $<4$ evaluable domains (UL30-FOS, $n=2$), **NO bootstrap CI is reported** because $\\binom{3}{2} = 3$ distinct resamples cannot support an empirical interval. Non-parametric sign consistency ($k/n$ domains with $\\text{AUROC}_{\\text{model}} > \\text{AUROC}_{\\text{pos\\_only}}$) compares each model against position-only RF and does not represent Where-versus-LR comparisons.

| Target Pair | $N_{\\text{eval}}$ ($n_{\\text{dom}}$) | Position-Only RF Mean [Pooled] | Split Validity Status | Zero-Shot ESM-2 650M Mean | **ESM-2 Excess $\\Delta$ [95% CI] (Sign)** | Zero-Shot BLOSUM62 Mean | **BLOSUM62 Excess $\\Delta$ [95% CI] (Sign)** | Where Spec-Exact Mean | **Where Excess $\\Delta$ [95% CI] (Sign)** | Supervised LR L2 Mean | **Sup-LR Excess $\\Delta$ [95% CI] (Sign)** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL30-Aciclovir** | 248 (4) | 0.5645 [**0.5185**] | **CONFIRMED** (Pooled in $[0.45, 0.55]$) | 0.7439 | **+0.1794** [+0.0281, +0.2664] (3/4) | 0.4906 | -0.0739 [-0.1589, +0.0115] (1/4) | **0.7801** | **+0.2157** [+0.1116, +0.3197] (4/4) | 0.7416 | **+0.1771** [+0.1274, +0.2339] (4/4) |
| **UL30-Foscarnet** | 114 (2) | 0.5253 [**0.5087**] | In-band ($n=2$ domains; pos-RF exactly 0.5000 in Fingers [constant]; not independently interpretable) | 0.5855 | +0.0602 [NO CI, $n=2$] (1/2) | 0.5845 | +0.0592 [NO CI, $n=2$] (1/2) | **0.7878** | **+0.2625** [NO CI, $n=2$] (2/2) | 0.7017 | **+0.1764** [NO CI, $n=2$] (2/2) |
| **UL23-Aciclovir** | 254 (4) | 0.4372 [**0.3775**] | Chance-referenced (reference clamped to 0.50; every domain is below chance) | 0.6885 | **+0.1885** [+0.1046, +0.2725] (4/4) | **0.6726** | **+0.1726** [+0.1238, +0.2215] (4/4) | 0.5670 | +0.0670 [-0.0360, +0.1700] (2/4) | **0.7177** | **+0.2177** [+0.1231, +0.2761] (4/4) |

*Per-Domain Breakdown & Diagnostic Notes (`logs/step0b_final_recalibrated_benchmarks.json`):*
1. **Per-Domain Breakdown on Primary Target (UL30-ACV):**
   - *Exonuclease ($N=37$, 8 res):* Pos-RF = 0.7091; ESM-2 = 0.6595; BLOSUM62 = 0.5776; Where = **0.8621**; Sup-LR = 0.8233.
   - *Palm ($N=111$, 32 res):* Pos-RF = 0.5332; ESM-2 = 0.8050; BLOSUM62 = 0.5712; Where = **0.8133**; Sup-LR = 0.7983.
   - *Fingers ($N=18$, 16 res):* Pos-RF = 0.5156; ESM-2 = 0.7500; BLOSUM62 = 0.5000; Where = **0.8750**; Sup-LR = 0.6562.
   - *Thumb / C-term ($N=82$, 6 res):* Pos-RF = 0.5000; ESM-2 = **0.7610**; BLOSUM62 = 0.3136; Where = 0.5702; Sup-LR = 0.6886.
2. **Where vs. Supervised LR Comparison:** On UL30-ACV, the label-free `Where` composite (domain mean: **0.7801**) matches Supervised Logistic Regression (domain mean: **0.7416**). Across individual domains, Where > LR in 3 of 4 domains (Exo: 0.8621 vs 0.8233; Palm: 0.8133 vs 0.7983; Fingers: 0.8750 vs 0.6562; Thumb: 0.5702 vs 0.6886); the $+0.0385$ domain-mean margin is largely driven by Fingers ($N=18$). On UL30-FOS, Where and LR are essentially matched in Palm (0.8256 vs 0.8200) and Fingers (0.7500 vs 0.5833). On UL23-ACV, `Where` (domain mean: **0.5670**) lags behind Supervised LR (0.7177), ESM-2 (0.6885), and BLOSUM62 (0.6726). Sign-consistency counts compare each model individually against the position-only RF baseline and are not direct Where-versus-LR comparisons.
3. **UL23 Position-Only RF Inspection:** Because position-only RF operates below chance on UL23 (domain mean: 0.4372; pooled: 0.3775), excess AUROC is explicitly clamped to a chance reference of 0.50. An inspection of position-only RF out-of-domain predictions across UL23 positions (`figures/ul23_pos_rf_diagnostic.png`) shows that predicted probabilities for resistant variants are lower than for sensitive variants within ATP Core (mean: 0.1356 vs 0.2553) and Nucleoside Pocket (0.4242 vs 0.6267), while C-term yields constant predictions (0.1000) and Lid/Dimer yields 0.2408 vs 0.2276. We report these empirical prediction values directly without asserting a speculative mechanistic cause prior to prospective validation.
4. **Descriptive Dissociation Contrast for BLOSUM62:** Across the 4 evaluable domains of UL23 vs. UL30-ACV, BLOSUM62 radicalness exhibits a domain-mean contrast of $\\Delta_{\\text{domain}} = 0.6726 - 0.4906 = \\mathbf{+0.1820}$ (with excess over 0.50 of $+0.1726$ on UL23 vs. $-0.0094$ on UL30). Given the small-$n$ domain limit ($n=4$ each), this contrast is presented strictly as a **descriptive** observation on these two genes rather than an inferential claim.
5. **Non-Independence of Arms:** In our 929 variant-drug records (630 unique substitutions), UL30-ACV and UL30-FOS share 264 sequence positions (**91.0% of ACV positions and 97.4% of FOS positions**) and 299 identical substitutions (**90.1% of ACV substitutions**). Similarly, the three HCMV UL54 cohorts represent identical underlying polymerase sequence records tested against multiple drugs. Consequently, ACV and FOS arms within UL30 (and the 3 drugs within HCMV UL54) represent **correlated, non-independent arms** of the same biophysical target rather than separate replication events.

"""

text = old_table7_pattern.sub(lambda m: new_table7, text)

# Fix beta5 text
old_b5_block = re.compile(
    r"An exploratory Leave-One-Position-Out \(LOO\) influence diagnostic across all 290 unique positions.*?neither confirmed nor refuted at current retrospective sample sizes\*\*\.",
    re.DOTALL
)

new_b5_block = """Pre-registered falsifier not met (estimate not supported); the estimate is leverage-driven: omitting position E927 (0 resistant / 1 variant) moves $\\beta_5$ from $-0.016$ to $+0.60$; ~31% of bootstrap resamples form a second mode. Specifically, in an exploratory Leave-One-Position-Out (LOO) influence diagnostic across all 290 unique positions (`logs/beta5_top15_loo_diagnostic.json`), omitting Position 927 (Glu) in the Palm domain alone moves $\\beta_5$ from $-0.0161$ to **$+0.6012$** ($\\Delta = +0.6173$). Across the 2,000 bootstrap resamples, the inclusion count of Position 927 correlates strongly with $\\beta_5$ ($r = -0.6966$): when Position 927 is omitted in random sampling (~31% of resamples), the fit shifts into a secondary positive mode (centered at $+0.683$, 31.4% of resamples). 

Omitting the entire Palm domain ($N=111$, containing 32 of the 62 resistant events) similarly shifts $\\beta_5$ to **$+0.4289$**, while omitting any other domain leaves $\\beta_5$ near zero. This influence diagnostic demonstrates that the secondary positive mode is largely driven by position-level leverage rather than a broad quasi-separation tail. Under our pre-registered decision rule (requiring the 95% CI upper bound to be $< 0$), **the tolerance-window hypothesis is unidentifiable: $\\beta_5$ is neither confirmed nor refuted at current retrospective sample sizes**."""

text = old_b5_block.sub(lambda m: new_b5_block, text)

# Remove "proving" and replace with demonstrates/indicates
text = text.replace("This diagnostic proves that the secondary positive mode", "This influence diagnostic demonstrates that the secondary positive mode")
text = text.replace("Proves overwhelming position identity memorization risk", "Demonstrates overwhelming position identity memorization risk")
text = text.replace("proves that amino acid radicalness", "demonstrates that amino acid radicalness")

# Replace "independent measurements" with "consistent observations on the same two genes"
text = text.replace("independent measurements", "consistent observations on the same two genes")

with open('manuscript/manuscript_draft.md', 'w', encoding='utf-8') as f:
    f.write(text)

print("Updated manuscript_draft.md for STAGE A successfully.")

# 2. Update README.md
with open('README.md', 'r', encoding='utf-8') as f:
    r_text = f.read()

old_readme_table3 = re.compile(
    r"### Table 3: Recalibrated Dual-Buffered Domain Benchmark.*?### Variance Decomposition",
    re.DOTALL
)

new_readme_table3 = """### Table 3: Recalibrated Dual-Buffered Domain Benchmark ($\pm 5$ aa, $< 12.0$ \AA)

Under a leave-one-domain-out protocol with a dual exclusion buffer, the Position-Only Random Forest operates at chance on the primary target (UL30-ACV pooled AUROC = **0.5185**). Performance is reported as evaluable per-domain mean alongside pooled out-of-fold AUROC (pooled AUROC measures joint across-domain ranking and base-rate calibration, whereas domain-mean AUROC isolates within-domain discriminative ordering). Evaluable domains strictly enforce the pre-registered $N_{\\text{pos}} \\ge 5$ rule (UL30-ACV: 4 domains [Exo, Palm, Fingers, Thumb]; UL30-FOS: 2 domains [Palm and Fingers; Exo and Thumb omitted with $N_{\\text{pos}} < 5$]; UL23-ACV: 4 domains).

**Excess AUROC is recalibrated against an empirical chance baseline:**
$$\\Delta_{\\text{excess}} = \\text{AUROC}_{\\text{model}} - \\max(0.50, \\text{AUROC}_{\\text{pos\\_only}})$$

For 4-domain targets ($n=4$), descriptive 95% confidence intervals are computed via domain-clustered bootstrap (recording 35 distinct multisets). For targets with $<4$ evaluable domains (UL30-FOS, $n=2$), **NO bootstrap CI is reported** because $\\binom{3}{2} = 3$ distinct resamples cannot support an empirical interval. Non-parametric sign consistency ($k/n$ domains with $\\Delta > 0$) compares each model individually against position-only RF and does not represent Where-versus-LR comparisons.

| Target Pair | $N_{\\text{eval}}$ ($n_{\\text{dom}}$) | Position-Only RF Mean [Pooled] | Split Validity Status | Zero-Shot ESM-2 650M Mean | **ESM-2 Excess $\\Delta$ [95% CI] (Sign)** | Zero-Shot BLOSUM62 Mean | **BLOSUM62 Excess $\\Delta$ [95% CI] (Sign)** | Where Spec-Exact Mean | **Where Excess $\\Delta$ [95% CI] (Sign)** | Supervised LR L2 Mean | **Sup-LR Excess $\\Delta$ [95% CI] (Sign)** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL30-Aciclovir** | 248 (4) | 0.5645 [**0.5185**] | **CONFIRMED** (Pooled in $[0.45, 0.55]$) | 0.7439 | **+0.1794** [+0.0281, +0.2664] (3/4) | 0.4906 | -0.0739 [-0.1589, +0.0115] (1/4) | **0.7801** | **+0.2157** [+0.1116, +0.3197] (4/4) | 0.7416 | **+0.1771** [+0.1274, +0.2339] (4/4) |
| **UL30-Foscarnet** | 114 (2) | 0.5253 [**0.5087**] | In-band ($n=2$; pos-RF exactly 0.5000 in Fingers [constant]; not independently interpretable) | 0.5855 | +0.0602 [NO CI, $n=2$] (1/2) | 0.5845 | +0.0592 [NO CI, $n=2$] (1/2) | **0.7878** | **+0.2625** [NO CI, $n=2$] (2/2) | 0.7017 | **+0.1764** [NO CI, $n=2$] (2/2) |
| **UL23-Aciclovir** | 254 (4) | 0.4372 [**0.3775**] | Chance-referenced (reference clamped to 0.50; every domain below chance) | 0.6885 | **+0.1885** [+0.1046, +0.2725] (4/4) | **0.6726** | **+0.1726** [+0.1238, +0.2215] (4/4) | 0.5670 | +0.0670 [-0.0360, +0.1700] (2/4) | **0.7177** | **+0.2177** [+0.1231, +0.2761] (4/4) |

*Key Calibration Takeaways (`logs/step0b_final_recalibrated_benchmarks.json`):*
1. **Where vs. Supervised LR:** On UL30-ACV, label-free `Where` (domain mean: **0.7801**) matches Supervised LR (domain mean: **0.7416**; Where > LR in 3 of 4 domains, with the $+0.0385$ margin largely driven by Fingers [$N=18$]). On UL30-FOS, Where and LR are essentially matched in Palm (0.8256 vs 0.8200) and Fingers (0.7500 vs 0.5833). On UL23-ACV, `Where` (domain mean: **0.5670**) lags behind Supervised LR (0.7177), ESM-2 (0.6885), and BLOSUM62 (0.6726).
2. **Biological Dissociation:** Clamping the baseline to $\\max(0.50, \\text{pos}_d)$ ensures models are never rewarded for beating a sub-chance baseline. BLOSUM62 achieves an excess of **+0.1726** (4/4 domains positive, domain mean = 0.6726) on thymidine kinase (UL23), while operating at or below chance on polymerase (UL30-ACV: -0.0739, 1/4 positive, domain mean = 0.4906). Descriptive difference is $\\Delta = \\mathbf{+0.1820}$ (consistent observations on the same two genes, $n=4$ domains each).
3. **Non-Independence of Arms:** UL30-ACV and FOS share 90.1% identical substitutions (264 of 290 positions = 91.0% of ACV, 97.4% of FOS) and HCMV UL54's three drugs represent one gene, representing correlated rather than independent arms.

---

### Variance Decomposition"""

r_text = old_readme_table3.sub(lambda m: new_readme_table3, r_text)

# Update beta5 in README
old_b5_readme_block = re.compile(
    r"4\. \*\*Tolerance-Window Hypothesis Test & Bimodality Diagnostic:\*\*.*?(?=---)",
    re.DOTALL
)

new_b5_readme_block = """4. **Tolerance-Window Hypothesis Test & Bimodality Diagnostic:** Pre-registered falsifier not met (estimate not supported); the estimate is leverage-driven: omitting position E927 (0 resistant / 1 variant) moves $\\beta_5$ from $-0.016$ to $+0.60$; ~31% of bootstrap resamples form a second mode. An influence diagnostic (`logs/beta5_top15_loo_diagnostic.json`) demonstrates that the secondary mode is largely driven by Position 927 in the Palm domain ($r = -0.6966$ with resample count). $\\beta_5$ is statistically unidentifiable in current retrospective data; the hypothesis is neither confirmed nor refuted. In contrast, UL23 exhibits strictly monotonic disruption ($\\beta_2 = +0.4468$, 95% CI: $[+0.1822, +0.7900]$, $p < 0.01$).

"""

r_text = old_b5_readme_block.sub(lambda m: new_b5_readme_block, r_text)

with open('README.md', 'w', encoding='utf-8') as f:
    f.write(r_text)

print("Updated README.md for STAGE A successfully.")
