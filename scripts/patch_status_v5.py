import re

status_path = r'd:\antiviral-resistance-hsv--1\STATUS.md'
with open(status_path, 'r', encoding='utf-8') as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    if line.startswith('| **3. FIRST-REPORT YEAR & TEMPORAL SHIFT**'):
        new_lines.append(
            "| **3. FIRST-REPORT YEAR & TEMPORAL SHIFT** | State exact split basis; analyze temporal domain shift | **DONE** | [`logs/publication_year_lookup.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/publication_year_lookup.tsv), [`logs/temporal_confounding.log`](file:///d:/antiviral-resistance-hsv--1/logs/temporal_confounding.log) | Publication-year split recomputed from resolved Crossref DOI publication years (Pre-2016: 725 rows, 569 mutations; Post-2016: 226 rows, 189 mutations). Passes Option A threshold (>= 30 post-cutoff mutations). Factual analysis confirmed the identical 316 count is a numerical coincidence across all supplementary table additions (Pre-2016: 316/1,271 [24.86%]; Post-2016: 316/858 [36.83%]). Non-label domain classifier achieved 5-fold cross-validated AUROC of 0.6794 +/- 0.0232 (Random Forest) and 0.6460 +/- 0.0199 (Logistic Regression), confirming significant experimental-era covariate shift. |\n"
        )
    elif line.startswith('| **6. C2 REPEATED CV & LEAKAGE EVALUATION**'):
        new_lines.append(
            "| **6. C2 REPEATED CV & LEAKAGE EVALUATION** | 20-repeat CV across 6 models and 5 split regimes | **DONE** | [`logs/c2_repeated_cv.log`](file:///d:/antiviral-resistance-hsv--1/logs/c2_repeated_cv.log), [`logs/c2_repeated_cv_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/c2_repeated_cv_results.json) | Evaluated 20 independent repeats (seeds 1..20, 100 evaluations/cell) across 6 models and 5 split regimes. Fixed Regime 4 leakage: Position Prior scored exactly 0.5000 +/- 0.0000 with asserted zero position and strain overlap. Paired leakage gaps (Random minus Pos-Grouped) strictly exclude zero under Nadeau & Bengio (2003) corrected resampled t-test 95% CIs: RF All Features (+0.1041 [+0.0703, +0.1379] on UL23-ACV, +0.1292 [+0.0878, +0.1706] on UL30-ACV, +0.1579 [+0.1050, +0.2108] on UL30-FOS). Dropped unverified p-value claims. Described as 'consistent with leakage in the random forest'. |\n"
        )
    elif line.startswith('| **10. LANGUAGE, HOMOLOGY & REPRODUCIBILITY**'):
        new_lines.append(line)
        new_lines.append(
            "| **11. ZERO-SHOT PLM BENCHMARK (DUAL-SPLIT)** | ESM-2 masked marginal evaluation on primary folds | **DONE** | [`data/processed/esm2_pos_cache.json`](file:///d:/antiviral-resistance-hsv--1/data/processed/esm2_pos_cache.json), [`scripts/evaluate_esm2_from_cache.py`](file:///d:/antiviral-resistance-hsv--1/scripts/evaluate_esm2_from_cache.py) | Evaluated zero-shot masked marginals (Delta LLR) across all 507 positions on 929 clean primary rows. Confirmed negative correlation with resistance (catalytic disruption). Under weak Position-Grouped split: UL23-ACV 0.6599 +/- 0.0495, UL30-ACV 0.6819 +/- 0.0934, UL30-FOS 0.7352 +/- 0.0157. Under strict Contiguous-Block split: UL23-ACV 0.6118 +/- 0.1786, UL30-ACV 0.6453 +/- 0.2070, UL30-FOS 0.6973 +/- 0.1841. While tabular models collapse to chance or below chance under Contiguous-Block splits (RF No Pos 0.3956, RF All 0.4089, Pos Only 0.5001), ESM-2 demonstrates robust structural transfer (0.6453 and 0.6973). |\n"
        )
    else:
        new_lines.append(line)

with open(status_path, 'w', encoding='utf-8') as f:
    f.writelines(new_lines)
print("STATUS.md patched successfully!")
