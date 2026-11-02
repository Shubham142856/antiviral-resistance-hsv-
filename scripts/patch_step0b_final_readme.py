"""
patch_step0b_final_readme.py
Updates README.md with final Step 0b benchmark numbers:
- Supervised LR column in Table 3
- NO bootstrap CI for FOS (n=2, 3 distinct resamples)
- 929 variant-drug records across 630 unique substitutions
- Non-independence of UL30-ACV/FOS (90.1% overlap) and HCMV UL54 (3 drugs, 1 gene)
- Beta5 LOO diagnostic (Pos 927 in Palm, r = -0.697)
"""
import re

with open('README.md', 'r', encoding='utf-8') as f:
    text = f.read()

# Replace substitutions count
text = text.replace(
    "929 single-site missense substitutions",
    "929 single-site missense variant-drug records (630 unique substitutions)"
)
text = text.replace(
    "929 single-site missense mutations",
    "929 single-site missense variant-drug records (630 unique amino acid substitutions)"
)

old_table3_pattern = re.compile(
    r"### Table 3: Recalibrated Dual-Buffered Domain Benchmark.*?### Variance Decomposition",
    re.DOTALL
)

new_table3 = """### Table 3: Recalibrated Dual-Buffered Domain Benchmark ($\pm 5$ aa, $< 12.0$ \AA)

Under a leave-one-domain-out protocol with a dual exclusion buffer, the Position-Only Random Forest operates at chance on the primary target (UL30-ACV pooled AUROC = **0.5185**). Evaluable domains strictly enforce the pre-registered $N_{\\text{pos}} \\ge 5$ rule (UL30-ACV: 4 domains [Exo, Palm, Fingers, Thumb]; UL30-FOS: 2 domains [Palm and Fingers; Exo and Thumb omitted with $N_{\\text{pos}} < 5$]; UL23-ACV: 4 domains [ATP Core, Nucleoside Pocket, Lid/Dimer, C-term]).

**Excess AUROC is recalibrated against an empirical chance baseline:**
$$\\Delta_{\\text{excess}} = \\text{AUROC}_{\\text{model}} - \\max(0.50, \\text{AUROC}_{\\text{pos\\_only}})$$

For 4-domain targets ($n=4$), descriptive 95% confidence intervals are computed via domain-clustered bootstrap (recording 35 distinct multisets). For targets with $<4$ evaluable domains (UL30-FOS, $n=2$), **NO bootstrap CI is reported** because $\\binom{3}{2} = 3$ distinct resamples cannot support an empirical interval. Non-parametric sign consistency ($k/n$ domains with $\\Delta > 0$) is reported for all targets.

| Target Pair | $N_{\\text{eval}}$ ($n_{\\text{dom}}$) | Position-Only RF Mean [Pooled] | Split Validity ($[0.45, 0.55]$) | Zero-Shot ESM-2 650M Mean | **ESM-2 Excess $\\Delta$ [95% CI] (Sign)** | Zero-Shot BLOSUM62 Mean | **BLOSUM62 Excess $\\Delta$ [95% CI] (Sign)** | Where Spec-Exact Mean | **Where Excess $\\Delta$ [95% CI] (Sign)** | Supervised LR L2 Mean | **Sup-LR Excess $\\Delta$ [95% CI] (Sign)** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL30-Aciclovir** | 248 (4) | 0.5645 [**0.5185**] | **CONFIRMED** (Pooled in $[0.45, 0.55]$) | 0.7439 | **+0.1794** [+0.0281, +0.2664] (3/4) | 0.4906 | -0.0739 [-0.1589, +0.0115] (1/4) | **0.7801** | **+0.2157** [+0.1116, +0.3197] (4/4) | 0.7416 | **+0.1771** [+0.1274, +0.2339] (4/4) |
| **UL30-Foscarnet** | 114 (2) | 0.5253 [**0.5087**] | **CONFIRMED** (In Palm & Fingers) | 0.5855 | +0.0602 [NO CI, $n=2$] (1/2) | 0.5845 | +0.0592 [NO CI, $n=2$] (1/2) | **0.7878** | **+0.2625** [NO CI, $n=2$] (2/2) | 0.7017 | **+0.1764** [NO CI, $n=2$] (2/2) |
| **UL23-Aciclovir** | 254 (4) | 0.4372 [**0.3775**] | **OUTSIDE** (Domain rate shift: 14%–75%) | 0.6885 | **+0.1885** [+0.1046, +0.2725] (4/4) | **0.6726** | **+0.1726** [+0.1238, +0.2215] (4/4) | 0.5670 | +0.0670 [-0.0360, +0.1700] (2/4) | **0.7177** | **+0.2177** [+0.1231, +0.2761] (4/4) |

*Key Calibration Takeaways (`logs/step0b_final_recalibrated_benchmarks.json`):*
1. **Where vs. Supervised LR:** On UL30-ACV, label-free `Where` (domain mean: **0.7801**) matches and marginally exceeds Supervised LR (domain mean: **0.7416**; Sign: 4/4 vs 4/4). However, on UL23-ACV, `Where` (domain mean: **0.5670**) lags behind Supervised LR (0.7177), ESM-2 (0.6885), and BLOSUM62 (0.6726), demonstrating that structural pocket prior is target-specific.
2. **Biological Dissociation:** Clamping the baseline to $\\max(0.50, \\text{pos}_d)$ ensures models are never rewarded for beating a sub-chance baseline. BLOSUM62 achieves an excess of **+0.1726** (4/4 domains positive, domain mean = 0.6726) on thymidine kinase (UL23), while operating at or below chance on polymerase (UL30-ACV: -0.0739, 1/4 positive, domain mean = 0.4906). Descriptive difference is $\\Delta = \\mathbf{+0.1820}$ (n=4 domains each).
3. **Non-Independence of Arms:** UL30-ACV and FOS share 90.1% identical substitutions (264 of 290 positions = 91.0% of ACV, 97.4% of FOS) and HCMV UL54's three drugs represent one gene, representing correlated rather than independent arms.

---

### Variance Decomposition"""

text = old_table3_pattern.sub(lambda m: new_table3, text)

# Update beta5 text in README
old_b5_readme = "4. **Tolerance-Window Hypothesis Test:** Pre-registered test of non-monotonic curvature on UL30 ($\beta_5$) yielded a null result ($\beta_5 = -0.0162$, 95% CI: $[-0.1876, +1.0402]$, spanning zero; not supported under pre-registered stopping rule). In contrast, UL23 exhibits strictly monotonic disruption ($\beta_2 = +0.4468$, 95% CI: $[+0.1822, +0.7900]$, $p < 0.01$)."

new_b5_readme = """4. **Tolerance-Window Hypothesis Test & Bimodality Diagnostic:** Pre-registered test of non-monotonic curvature on UL30 ($\beta_5$) yielded an unidentifiable result ($\beta_5 = -0.0162$, 95% CI: $[-0.1876, +1.0402]$). An exploratory Leave-One-Position-Out diagnostic (`logs/beta5_top15_loo_diagnostic.json`) demonstrates that the bimodal bootstrap distribution (Mode 1: 68.5% at $-0.018$; Mode 2: 31.4% at $+0.683$) is driven by high leverage at **Position 927 (Glu)** in the Palm domain (LOO $\beta_5 = +0.6012$, $r = -0.6966$ with resample count). $\beta_5$ is statistically unidentifiable in current retrospective data; the hypothesis is neither confirmed nor refuted. In contrast, UL23 exhibits strictly monotonic disruption ($\beta_2 = +0.4468$, 95% CI: $[+0.1822, +0.7900]$, $p < 0.01$)."""

text = text.replace(old_b5_readme, new_b5_readme)

with open('README.md', 'w', encoding='utf-8') as f:
    f.write(text)

print("Updated README.md successfully.")
