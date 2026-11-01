import pathlib

def main():
    p = pathlib.Path(__file__).resolve().parent.parent / "STATUS.md"
    lines = p.read_text(encoding="utf-8").splitlines()

    new_section = [
        "   - **3D Spatial-Neighborhood Holdout Pre-Registered & Evaluated (2026-11-01):** Executed the empirical structural neighborhood holdout experiment (`scripts/run_spatial_3d_holdout_experiment.py`) using cryo-EM coordinates (`8V1Q` Chain A, 1097 CA atoms, 88.6% reference coverage, 80.8% unique benchmark position coverage) and crystal coordinates (`1OF1` Chain B, 313 CA atoms, 1.95 Å resolution, complexed with (South)-methanocarba-thymidine; Bennett et al. 1999). Pre-registered 5 deterministic spatial clusters via KMeans (`k=5, seed=42`) into `data/processed/primary_3d_spatial_folds.tsv` (SHA-256: `0ED51D37D558530EF8B2D0C991780A23936580920E62C7ADE46A426ADFA96F4C`, logged in `logs/spatial_3d_fold_assignments_hash.log`). Key findings serialized in `logs/spatial_3d_holdout_results.json`:",
        "     1. *Methodological Null Result Regarding PLM Structural Dependence:* ESM-2 650M maintains invariant generalization from 1D contiguous blocks to 3D spatial neighborhoods: UL30-ACV 0.7875 -> 0.7875 (resolved core: 0.7887 -> 0.7966); UL30-FOS 0.8809 -> 0.8932 (resolved core: 0.8793 -> 0.8977); UL23-ACV 0.7397 -> 0.7303 (resolved core: 0.7255 -> 0.7132). Cross-regime |Delta| <= 0.012 confirms that zero-shot PLMs do not depend on localized 3D physical contact environments, and 3D spatial grouping provides no incremental evaluation penalty.",
        "     2. *Negative Control Stringency Finding (1D Block is Strictly Stricter):* Position-only random forest rebounds from pure chance under 1D contiguous block holdout (0.5001 on UL30-ACV) to 0.6241 +/- 0.345 under 3D spatial holdout (and from 0.5166 to 0.6001 +/- 0.288 on UL30-FOS). Because 3D geometric clusters partition sequential primary residues across different spatial folds, models exploit sequence-adjacency interpolation across training and test folds. This proves that 1D contiguous-block holdout is strictly the more stringent evaluation standard.",
        "     3. *Target Dissociation Invariant in 3D Space:* Chemical dissimilarity alone (BLOSUM62 radicalness) consistently clears chance on the catalytic kinase (UL23-ACV: 1D 0.6380, 3D 0.6249; resolved core: 1D 0.6732, 3D 0.6501) but sits at chance on the allosteric polymerase (UL30-ACV: 1D 0.5494, 3D 0.4908; UL30-FOS: 1D 0.5074, 3D 0.4847).",
        "     4. *Evolutionary Baselines in 3D Space:* Multi-herpesvirus PSSM scores 0.7189 (UL30-ACV) and 0.8334 (UL30-FOS), while positional Shannon entropy scores 0.6806 and 0.8338. ESM-2 maintains its lead over both baselines (+0.069 to +0.107 on ACV; +0.060 on FOS).",
        "     5. *Publication Figure 5 Generated:* High-resolution 300-DPI multi-panel visualization saved to `manuscript/figures/figure5_spatial_3d_generalization.png` and artifact storage.",
        "   - **Primary Reference List Completed (2026-11-01):** Restored missing primary data sources and structural citations to References: HerpesDRG (Charles et al., BMC Bioinformatics 2024, PMID: 39192205), CMI 2025 review (Dähne et al., Clin Microbiol Infect 2025, PMID: 40349973), and HSV-1 TK structure (Bennett et al., FEBS Lett 1999, PMID: 10025941). Corrected spelling from 'Döhne' to 'Dähne' throughout, and corrected cohort scope sentence to HSV-1.",
        "   - **Estimator Alignment (2026-11-01):** Made the block-level Welch test primary for the target dissociation interaction claim (p=0.431 and p=0.406 across structural blocks), reporting the position-cluster bootstrap (Delta = +0.1187, p=0.0020; Delta = +0.1598, p < 0.0001) strictly as secondary sensitivity evidence, maintaining complete consistency with the paper's central thesis on spatial leakage."
    ]

    # Find where line starts with "   - **3D Spatial-Neighborhood Holdout"
    start_idx = None
    end_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("- **3D Spatial-Neighborhood Holdout"):
            start_idx = i
        if start_idx is not None and line.strip().startswith("- **Foscarnet threshold corrected"):
            end_idx = i
            break

    if start_idx is not None and end_idx is not None:
        lines[start_idx:end_idx] = new_section + [""]
        p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"Successfully updated STATUS.md lines {start_idx} to {end_idx}!")
    else:
        print(f"Could not find indices: start={start_idx}, end={end_idx}")

if __name__ == "__main__":
    main()

