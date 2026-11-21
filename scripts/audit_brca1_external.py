"""Apply the prespecified position-only audit to the Findlay BRCA1 assay."""

import hashlib
import json
import re
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

from audit_position_control import evaluate, paired_position_interval
from matplotlib import pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SCORES = ROOT / "data/brca1_findlay_scores.csv"
TABLE = ROOT / "data/GSE117159_Supplementary_Table_1.xlsx"
PLAN = ROOT / "docs/BRCA1_EXTERNAL_AUDIT.md"


def annotations():
    """Read the first worksheet with its row-four column names and shared strings."""
    with ZipFile(TABLE) as book:
        strings = ["".join(item.itertext()) for item in
                   ElementTree.fromstring(book.read("xl/sharedStrings.xml"))]
        sheet = ElementTree.fromstring(book.read("xl/worksheets/sheet1.xml"))
    records = []
    for row in sheet.findall(".//{*}row"):
        number = int(row.attrib["r"])
        if number < 4:
            continue
        cells = {}
        for cell in row:
            key = re.match(r"[A-Z]+", cell.attrib["r"]).group()
            value = cell.findtext("{*}v")
            if value is not None:
                cells[key] = strings[int(value)] if cell.get("t") == "s" else value
        if number == 4:
            headers = cells
        else:
            records.append({headers[key]: value for key, value in cells.items()})
    return pd.DataFrame(records)


def adjacent_site_test(frame):
    sites = frame.groupby("pos").agg(rate=("target_binary", "mean"),
                                     exon=("experiment", "first"),
                                     exon_count=("experiment", "nunique"))
    groups = []
    by_exon = {}
    for exon, group in sites[sites.exon_count == 1].groupby("exon"):
        edges = np.flatnonzero(np.diff(group.index) == 1)
        if len(edges):
            rates = group.rate.to_numpy()
            groups.append((rates, edges))
            # Every edge draws two distinct sites uniformly under the shuffle.
            i, j = np.triu_indices(len(rates), 1)
            by_exon[exon] = dict(pairs=len(edges),
                                 observed=np.abs(rates[edges + 1] - rates[edges]).mean(),
                                 expected=np.abs(rates[i] - rates[j]).mean())
    pairs = sum(len(edges) for _, edges in groups)
    rng = np.random.default_rng(42)

    def difference(shuffle):
        total = 0.
        for rates, edges in groups:
            values = rng.permutation(rates) if shuffle else rates
            total += np.abs(values[edges + 1] - values[edges]).sum()
        return total / pairs

    observed = difference(False)
    null = np.array([difference(True) for _ in range(9999)])
    excess = sum(e["pairs"] * (e["expected"] - e["observed"]) for e in by_exon.values())
    leave_one_out = {key: (excess - e["pairs"] * (e["expected"] - e["observed"])) /
                     (pairs - e["pairs"]) for key, e in by_exon.items()}
    return dict(pairs=pairs, exons=len(groups), ambiguous_sites=int((sites.exon_count > 1).sum()),
                observed=observed, null_mean=null.mean(), permutations=len(null),
                per_exon=by_exon, leave_one_exon_out_expected_minus_observed=leave_one_out,
                lower_tail_p=(1 + np.count_nonzero(null <= observed + 1e-12)) / (1 + len(null)))


def plot_exons(result):
    exons = sorted(result["per_exon"], key=lambda name: int(name.lstrip("X")))
    fig, ax = plt.subplots(figsize=(7.6, 3.7))
    observed = [result["per_exon"][name]["observed"] for name in exons]
    expected = [result["per_exon"][name]["expected"] for name in exons]
    x = np.arange(len(exons))
    ax.vlines(x, observed, expected, color="#bbc5ca")
    ax.scatter(x, expected, facecolors="white", edgecolors="#296778", label="Within-exon expectation")
    ax.scatter(x, observed, marker="D", color="#ab442b", label="Adjacent sites")
    ax.set(xticks=x, xticklabels=[name.lstrip("X") for name in exons],
           xlabel="BRCA1 exon experiment", ylabel="Mean absolute difference in\nnonfunctional fraction", ylim=(0, None))
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.22))
    fig.tight_layout()
    fig.savefig(ROOT / "manuscript/figures/brca1_exon_adjacency.png", dpi=300)
    plt.close(fig)


def main():
    scores = pd.read_csv(SCORES)
    meta = annotations()
    meta["hgvs_nt"] = meta.transcript_ID + ":" + meta.transcript_variant
    joined = scores.merge(meta, on="hgvs_nt", how="outer", validate="one_to_one", indicator=True)
    if not joined["_merge"].eq("both").all():
        raise ValueError("GEO and MaveDB disagree on variant membership")
    error = (joined.score - joined["function.score.mean"].astype(float)).abs().max()
    if error > 1e-8:
        raise ValueError(f"Function scores disagree by up to {error}")
    missense = joined[joined.consequence.str.casefold().eq("missense")]
    frame = missense[missense["func.class"].isin(["FUNC", "LOF"])].copy()
    frame["pos"] = frame.aa_pos.astype(int)
    frame = frame.sort_values(["pos", "hgvs_nt"]).reset_index(drop=True)
    frame["target_binary"] = frame["func.class"].eq("LOF").astype(int)
    frame["norm_pos"] = frame.pos / 1863
    for fold, (_, test) in enumerate(GroupKFold(5).split(frame, groups=frame.pos)):
        frame.loc[test, "grouped"] = fold
    for fold, sites in enumerate(np.array_split(np.unique(frame.pos), 5)):
        frame.loc[frame.pos.isin(sites), "blocked"] = fold
    grouped, group_scores = evaluate(frame, "grouped")
    blocked, block_scores = evaluate(frame, "blocked")
    out = dict(hashes={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (SCORES, TABLE, PLAN)}, joined=len(joined), max_score_difference=error,
               missense_classes=missense["func.class"].value_counts().to_dict(),
               grouped=grouped, blocked=blocked, paired=paired_position_interval(frame, group_scores, block_scores),
               adjacency=adjacent_site_test(frame))
    path = ROOT / "logs/brca1_external_audit.json"
    path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    plot_exons(out["adjacency"])
    print(json.dumps({key: out[key] for key in ("paired", "adjacency")}, indent=2))


if __name__ == "__main__":
    main()
