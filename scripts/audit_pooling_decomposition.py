"""Decompose pooled out-of-fold AUROC into within-fold and between-fold pair concordance.

Reads only existing audit logs (no refitting):
  logs/position_control_folds.json   - HSV position-only control, 3 arms x {grouped, 3D, block}
  logs/brca1_external_audit.json     - BRCA1 position-only control, {grouped, block}

For every target/cohort/design it recovers the identity

    pooled = w_within * C_within + w_between * C_between ,
    w_within = n_within / (n_within + n_between)

verifies it against the logged pooled AUROC, and reports

    offset = pooled - mean_within_fold_AUROC          (pooling offset)
    gap    = w_between * (C_between - C_within)        (between-pair contribution)

It also decomposes the block -> 3D change in pooled AUROC into the part carried by
within-fold pairs and the part carried by between-fold pairs, using the Oaxaca-style
split evaluated at the block-design weights:

    d_pooled = d_within_component + d_between_component + d_weight_component

Writes logs/pooling_decomposition.json and manuscript/figures/pooling_decomposition.png.
"""

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
LOGS = ROOT / "logs"
FIGS = ROOT / "manuscript" / "figures"

TOL = 5e-4  # identity tolerance in AUROC units


def decompose(entry):
    """Return the within/between decomposition of one logged design."""
    ps = entry["pair_summary"]
    n_w = ps["within"]["count"]
    n_b = ps["between"]["count"]
    c_w = ps["within"]["concordance"]
    c_b = ps["between"]["concordance"]
    n = n_w + n_b
    w_w = n_w / n
    w_b = n_b / n
    reconstructed = w_w * c_w + w_b * c_b
    pooled = entry["pooled"]
    folds = [f for f in entry["folds"] if f.get("auc") is not None]
    return {
        "records": entry["records"],
        "positions": entry["positions"],
        "median_gap": entry["median_gap"],
        "pairs_within": n_w,
        "pairs_between": n_b,
        "frac_between": w_b,
        "concordance_within": c_w,
        "concordance_between": c_b,
        "pooled_logged": pooled,
        "pooled_reconstructed": reconstructed,
        "identity_abs_error": abs(reconstructed - pooled),
        "identity_holds": abs(reconstructed - pooled) < TOL,
        "mean_within_fold_auroc": entry["within_fold_mean"],
        "evaluable_folds": len(folds),
        "total_folds": len(entry["folds"]),
        "pooling_offset_vs_fold_mean": pooled - entry["within_fold_mean"],
        "between_pair_contribution": w_b * (c_b - c_w),
        "prevalence_only_pooled": entry.get("prevalence_only_pooled"),
    }


def change_decomposition(base, alt):
    """Split pooled(alt) - pooled(base) into within, between and weight components."""
    wb_base = base["frac_between"]
    ww_base = 1.0 - wb_base
    d_within = ww_base * (alt["concordance_within"] - base["concordance_within"])
    # components evaluated at the base design's pair weights
    d_between = wb_base * (alt["concordance_between"] - base["concordance_between"])
    d_total = alt["pooled_logged"] - base["pooled_logged"]
    d_weight = d_total - d_within - d_between
    return {
        "d_pooled": d_total,
        "d_within_component": d_within,
        "d_between_component": d_between,
        "d_weight_component": d_weight,
        "between_share_of_change": (d_between / d_total) if d_total != 0 else None,
        "d_mean_within_fold_auroc": alt["mean_within_fold_auroc"] - base["mean_within_fold_auroc"],
    }


def main():
    hsv = json.loads((LOGS / "position_control_folds.json").read_text())
    brca = json.loads((LOGS / "brca1_external_audit.json").read_text())

    out = {
        "identity": "pooled_AUROC = w_within * C_within + w_between * C_between",
        "identity_tolerance": TOL,
        "sources": {
            "hsv": "logs/position_control_folds.json",
            "hsv_source_sha256": hsv["source_sha256"],
            "sklearn": hsv["sklearn"],
            "brca1": "logs/brca1_external_audit.json",
        },
        "designs": {},
        "changes": {},
    }

    for target, cohorts in hsv["targets"].items():
        for cohort, designs in cohorts.items():
            for design, entry in designs.items():
                key = f"HSV|{target}|{cohort}|{design}"
                out["designs"][key] = decompose(entry)

    for design, entry in brca.items():
        if design not in ("grouped", "blocked"):
            continue
        name = "Position-Grouped" if design == "grouped" else "Contiguous Block"
        out["designs"][f"BRCA1|BRCA1-SGE|all|{name}"] = decompose(entry)

    # Block -> 3D change, coordinate-matched cohort (the comparison the project
    # previously interpreted as spatial-adjacency leakage).
    for target in hsv["targets"]:
        cm = hsv["targets"][target]["coordinate_matched"]
        base = out["designs"][f"HSV|{target}|coordinate_matched|Contiguous Block"]
        alt = out["designs"][f"HSV|{target}|coordinate_matched|3D k-Means Spatial"]
        out["changes"][f"{target}|block_to_3D"] = change_decomposition(base, alt)
        grp = out["designs"][f"HSV|{target}|coordinate_matched|Position-Grouped"]
        out["changes"][f"{target}|block_to_grouped"] = change_decomposition(base, grp)
        assert cm is not None

    # How much of the block -> 3D pooled change is reproduced by the
    # training-prevalence-only predictor? That predictor is constant inside each
    # test fold, so every within-fold pair ties and every evaluable fold has
    # AUROC exactly 0.5 by construction. Any movement in its pooled score is
    # carried entirely by between-fold pairs.
    out["prevalence_only_changes"] = {}
    for target in hsv["targets"]:
        d = out["designs"]
        blk = d[f"HSV|{target}|coordinate_matched|Contiguous Block"]
        spa = d[f"HSV|{target}|coordinate_matched|3D k-Means Spatial"]
        prev_change = spa["prevalence_only_pooled"] - blk["prevalence_only_pooled"]
        fitted_change = spa["pooled_logged"] - blk["pooled_logged"]
        out["prevalence_only_changes"][f"{target}|block_to_3D"] = {
            "prevalence_only_block": blk["prevalence_only_pooled"],
            "prevalence_only_3D": spa["prevalence_only_pooled"],
            "prevalence_only_change": prev_change,
            "fitted_position_only_change": fitted_change,
            "share_of_fitted_change": prev_change / fitted_change if fitted_change else None,
        }

    failures = [k for k, v in out["designs"].items() if not v["identity_holds"]]
    out["identity_failures"] = failures

    (LOGS / "pooling_decomposition.json").write_text(json.dumps(out, indent=2))

    # ---- report ----
    print(f"{'design':52s} {'n_btw%':>7s} {'C_w':>7s} {'C_b':>7s} {'pooled':>7s} {'recon':>7s} {'foldmu':>7s} {'offset':>8s}")
    for k, v in out["designs"].items():
        print(
            f"{k:52s} {100*v['frac_between']:7.1f} {v['concordance_within']:7.3f} "
            f"{v['concordance_between']:7.3f} {v['pooled_logged']:7.3f} "
            f"{v['pooled_reconstructed']:7.3f} {v['mean_within_fold_auroc']:7.3f} "
            f"{v['pooling_offset_vs_fold_mean']:+8.3f}"
        )
    print(f"\nidentity failures: {failures if failures else 'none'}")
    print("\nblock -> alternative design, decomposition of the pooled change:")
    for k, v in out["changes"].items():
        share = v["between_share_of_change"]
        share_s = f"{100*share:.0f}%" if share is not None else "n/a"
        print(
            f"  {k:40s} d_pooled={v['d_pooled']:+.3f}  within={v['d_within_component']:+.3f}  "
            f"between={v['d_between_component']:+.3f}  weight={v['d_weight_component']:+.3f}  "
            f"between share={share_s}  d_foldmean={v['d_mean_within_fold_auroc']:+.3f}"
        )

    print("\ntraining-prevalence-only predictor (constant within each fold, fold AUROC = 0.5 exactly):")
    for k, v in out["prevalence_only_changes"].items():
        print(
            f"  {k:30s} block={v['prevalence_only_block']:.3f} -> 3D={v['prevalence_only_3D']:.3f} "
            f"({v['prevalence_only_change']:+.3f}); fitted {v['fitted_position_only_change']:+.3f}; "
            f"share={100*v['share_of_fitted_change']:.0f}%"
        )

    # ---- figure: pooled vs mean within-fold AUROC, by design ----
    order = [
        ("UL30-Aciclovir", "UL30 / Aciclovir"),
        ("UL30-Foscarnet", "UL30 / Foscarnet"),
        ("UL23-Aciclovir", "UL23 / Aciclovir"),
    ]
    designs = ["Contiguous Block", "3D k-Means Spatial", "Position-Grouped"]
    labels = ["Contiguous\nblock", "3D spatial\nclusters", "Position\ngrouped"]
    fig, axes = plt.subplots(1, 4, figsize=(13.5, 3.9), sharey=True)
    for ax, (target, title) in zip(axes, order):
        pooled = [out["designs"][f"HSV|{target}|coordinate_matched|{d}"]["pooled_logged"] for d in designs]
        foldmu = [
            out["designs"][f"HSV|{target}|coordinate_matched|{d}"]["mean_within_fold_auroc"] for d in designs
        ]
        x = range(len(designs))
        ax.axhline(0.5, color="0.6", ls=":", lw=1)
        ax.plot(x, pooled, "D-", color="#b2182b", ms=8, label="Pooled out-of-fold")
        ax.plot(x, foldmu, "o-", color="#2166ac", ms=8, label="Mean within-fold")
        for xi, (p, f) in enumerate(zip(pooled, foldmu)):
            ax.annotate("", xy=(xi, p), xytext=(xi, f), arrowprops=dict(arrowstyle="-", color="0.75", lw=1))
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.set_ylim(0.35, 0.95)
    ax = axes[3]
    bd = ["Contiguous Block", "Position-Grouped"]
    pooled = [out["designs"][f"BRCA1|BRCA1-SGE|all|{d}"]["pooled_logged"] for d in bd]
    foldmu = [out["designs"][f"BRCA1|BRCA1-SGE|all|{d}"]["mean_within_fold_auroc"] for d in bd]
    ax.axhline(0.5, color="0.6", ls=":", lw=1)
    ax.plot([0, 2], pooled, "D-", color="#b2182b", ms=8)
    ax.plot([0, 2], foldmu, "o-", color="#2166ac", ms=8)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Contiguous\ngroups", "", "Position\ngrouped"], fontsize=8)
    ax.set_title("BRCA1 SGE assay", fontsize=10)
    ax.set_xlim(-0.5, 2.5)
    axes[0].set_ylabel("AUROC, position-only control")
    axes[0].legend(fontsize=8, loc="upper left", frameon=False)
    fig.suptitle(
        "Pooled and within-fold estimates of the same predictions diverge in a partition-dependent direction",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "pooling_decomposition.png", dpi=200)
    print(f"\nwrote {LOGS/'pooling_decomposition.json'}")
    print(f"wrote {FIGS/'pooling_decomposition.png'}")


if __name__ == "__main__":
    main()
