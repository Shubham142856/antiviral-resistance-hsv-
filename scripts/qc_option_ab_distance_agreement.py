"""
qc_option_ab_distance_agreement.py
==================================
STEP 1 label-free diagnostic for the Where-components ablation.

Question: is the Option A -> Option B distance-only AUROC swing a coordinate-provenance
effect, or an artefact of the reference atom set / unmapped-position handling?

As executed in scripts/run_step2_hcmv_replication.py the two options differ in TWO ways
at once:
  Option A : Swiss-Model UL54 coords   + all atoms of the 19-residue HCMV active-site shell
  Option B : 8EXX HSV-1 coords         + PPF/MG ligand atoms + residues {717,888,889} only
This script decomposes those two changes by harmonizing the reference set across coordinate
sources, and adds an independent experimental coordinate replicate (8V1T).

CRITICAL CONSTRAINT: strictly zero resistance labels, fold-changes or phenotypes are read.
Only variant POSITIONS are taken from HerpesDRG, to define the evaluated position set.
"""

import json
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser, PDBParser
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

HSV1_SHELL_19 = [717, 718, 719, 720, 721, 722, 723, 785, 789, 811,
                 812, 815, 816, 818, 819, 887, 888, 925, 927]
HSV1_TRIAD = [717, 888, 889]


def min_dist(ref_atoms, ca_coords):
    ref = np.asarray(ref_atoms)
    return {p: float(np.min(np.linalg.norm(ref - c, axis=1))) for p, c in ca_coords.items()}


def collect(model, chain_id, residues=(), ligands=()):
    atoms = []
    for res in model.get_residues():
        if res.get_resname() in ligands:
            atoms += [a.get_coord() for a in res]
        if res.get_parent().id == chain_id and res.id[0] == " " and res.id[1] in residues:
            atoms += [a.get_coord() for a in res]
    return atoms


def ca_of(model, chain_id=None):
    out = {}
    for chain in model:
        if chain_id is not None and chain.id != chain_id:
            continue
        for res in chain:
            if "CA" in res:
                out[res.id[1]] = res["CA"].get_coord()
    return out


def main():
    sys.path.insert(0, ROOT)
    from scripts.map_active_site_to_hcmv import get_gene_seq

    print("=" * 78)
    print("STEP 1  LABEL-FREE: Option A vs Option B active-site distance agreement")
    print("=" * 78)

    # ---- alignment HSV-1 UL30 <-> HCMV UL54 (same settings as the runner) ----
    seq_hsv = get_gene_seq("NC_001806.2", "UL30")
    seq_hcmv = get_gene_seq("NC_006273.2", "UL54")
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -10.0, -0.5, "global"
    aln = al.align(seq_hsv, seq_hcmv)[0]
    hi = ci = 0
    hcmv_to_hsv, hsv_to_hcmv = {}, {}
    for r, q in zip(aln[0], aln[1]):
        if r != "-":
            hi += 1
        if q != "-":
            ci += 1
        if r != "-" and q != "-":
            hcmv_to_hsv[ci] = hi
            hsv_to_hcmv[hi] = ci
    print("")
    print("Alignment: HSV-1 UL30 %d aa <-> HCMV UL54 %d aa; %d aligned pairs"
          % (len(seq_hsv), len(seq_hcmv), len(hcmv_to_hsv)))

    # ---- coordinate sources ----
    sm = list(PDBParser(QUIET=True).get_structure(
        "UL54", os.path.join(ROOT, "data", "structures",
                             "P08546_swissmodel_8oj7.pdb")).get_models())[0]
    ca_sm = ca_of(sm)

    cif = MMCIFParser(QUIET=True)
    exx = list(cif.get_structure(
        "8EXX", os.path.join(ROOT, "data", "structures", "8exx.cif")).get_models())[0]
    v1t = list(cif.get_structure(
        "8V1T", os.path.join(ROOT, "data", "structures", "8v1t.cif")).get_models())[0]
    ca_exx, ca_v1t = ca_of(exx, "A"), ca_of(v1t, "A")

    with open(os.path.join(ROOT, "logs", "hcmv_ul54_active_site_mapping.json")) as f:
        hcmv_shell_19 = sorted(set(json.load(f)["hcmv_active_site_positions"]))
    hcmv_triad = sorted({hsv_to_hcmv[p] for p in HSV1_TRIAD if p in hsv_to_hcmv})
    print("HCMV 19-residue shell : %s" % hcmv_shell_19)
    print("HCMV triad (HSV-1 %s projected) : %s" % (HSV1_TRIAD, hcmv_triad))

    # ---- variants: all keyed in HCMV UL54 numbering ----
    variants = {}
    sm_chain = list(sm)[0].id
    variants["A_as_executed (SwissModel + shell19)"] = min_dist(
        collect(sm, sm_chain, residues=hcmv_shell_19), ca_sm)
    variants["A_triad (SwissModel + triad)"] = min_dist(
        collect(sm, sm_chain, residues=hcmv_triad), ca_sm)

    for tag, mdl, ca, lig, src in [("B", exx, ca_exx, ("PPF", "MG"), "8EXX"),
                                   ("C", v1t, ca_v1t, ("AVP", "MG"), "8V1T")]:
        d_lt = min_dist(collect(mdl, "A", residues=HSV1_TRIAD, ligands=lig), ca)
        d_sh = min_dist(collect(mdl, "A", residues=HSV1_SHELL_19), ca)
        name_lt = ("B_as_executed (8EXX + ligand+triad)" if tag == "B"
                   else "C_ligand_triad (8V1T)")
        variants[name_lt] = {p: d_lt[h] for p, h in hcmv_to_hsv.items() if h in d_lt}
        variants["%s_shell19 (%s + shell19)" % (tag, src)] = {
            p: d_sh[h] for p, h in hcmv_to_hsv.items() if h in d_sh}

    # ---- evaluated position set: LABEL-FREE (positions only) ----
    raw = pd.read_csv(os.path.join(ROOT, "data", "raw", "herpesdrg-db", "herpesdrg-db.tsv"),
                      sep="\t", usecols=["virus", "gene", "aa_change"])
    ul54 = raw[(raw.virus == "HCMV") & (raw.gene == "UL54")].copy()
    pos = [int(m.group(2)) for m in (re.match(r"^([A-Z])(\d+)([A-Z])$", str(s).strip())
                                     for s in ul54.aa_change) if m]
    cohort = sorted(set(pos))
    print("")
    print("LABEL-FREE cohort: %d UL54 rows -> %d parsed missense -> %d unique positions"
          % (len(ul54), len(pos), len(cohort)))

    # ---- coverage / unmapped handling ----
    print("")
    print("-" * 78)
    print("COVERAGE over evaluated positions (runner imputes missing with max distance)")
    print("-" * 78)
    print("%-46s %9s %8s %9s" % ("variant", "covered", "missing", "impute"))
    miss = {}
    for k, v in variants.items():
        have = [p for p in cohort if p in v and not np.isnan(v[p])]
        gone = [p for p in cohort if p not in have]
        miss[k] = gone
        mx = np.nanmax(list(v.values()))
        print("%-46s %9d %8d %9.2f" % (k, len(have), len(gone), mx))
    for k, g in miss.items():
        if g:
            print("")
            print("  MISSING in %s (n=%d): %s" % (k, len(g), g))

    # ---- pairwise agreement on positions covered by BOTH (no imputation) ----
    print("")
    print("-" * 78)
    print("PAIRWISE AGREEMENT (common covered positions, imputed values excluded)")
    print("-" * 78)
    print("%-70s %4s %9s %8s %8s" % ("pair", "n", "spearman", "med|d|", "max|d|"))
    keys = list(variants)
    pairs_out = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = variants[keys[i]], variants[keys[j]]
            common = [p for p in cohort
                      if p in a and p in b and not np.isnan(a[p]) and not np.isnan(b[p])]
            if len(common) < 5:
                continue
            x = np.array([a[p] for p in common])
            y = np.array([b[p] for p in common])
            rs = spearmanr(x, y).statistic
            d = np.abs(x - y)
            print("%-34s vs %-34s %4d %9.3f %8.2f %8.2f"
                  % (keys[i][:33], keys[j][:33], len(common), rs, np.median(d), np.max(d)))
            pairs_out.append({"a": keys[i], "b": keys[j], "n_common": len(common),
                              "spearman": round(float(rs), 4),
                              "median_abs_delta_angstrom": round(float(np.median(d)), 3),
                              "max_abs_delta_angstrom": round(float(np.max(d)), 3)})

    # ---- ligand reference-set geometric coherence (label-free) ----
    lig_audit = {}
    for src, mdl, ligs in [("8EXX", exx, ("PPF", "MG")), ("8V1T", v1t, ("AVP", "MG"))]:
        pts, members = [], []
        for res in mdl.get_residues():
            if res.get_resname() in ligs:
                members.append({"resname": res.get_resname(), "chain": res.get_parent().id,
                                "resseq": res.id[1],
                                "centroid": [round(float(z), 1)
                                             for z in np.mean([a.get_coord() for a in res], axis=0)]})
                pts += [a.get_coord() for a in res]
        pts = np.asarray(pts)
        tri = np.asarray(collect(mdl, "A", residues=HSV1_TRIAD))
        lig_audit[src] = {
            "ligand_members": members,
            "n_atoms": int(len(pts)),
            "max_pairwise_spread_angstrom": round(float(
                np.max(np.linalg.norm(pts[:, None] - pts[None], axis=-1))), 1),
            "max_ligand_atom_to_chainA_triad_centroid_angstrom": round(float(
                np.max(np.linalg.norm(pts - tri.mean(0), axis=1))), 1)}

    out = {
        "scope": "LABEL-FREE: no resistance labels, fold-changes or phenotypes read",
        "evaluated_positions": {"source": "data/raw/herpesdrg-db/herpesdrg-db.tsv",
                                "filter": "virus==HCMV & gene==UL54, single missense",
                                "n_rows": int(len(ul54)), "n_parsed": len(pos),
                                "n_unique_positions": len(cohort)},
        "alignment": {"hsv1_ul30_aa": len(seq_hsv), "hcmv_ul54_aa": len(seq_hcmv),
                      "n_aligned_pairs": len(hcmv_to_hsv)},
        "reference_sets": {"hcmv_shell_19": hcmv_shell_19, "hcmv_triad": hcmv_triad,
                           "hsv1_shell_19": HSV1_SHELL_19, "hsv1_triad": HSV1_TRIAD},
        "coverage": {k: {"covered": int(len([p for p in cohort
                                             if p in v and not np.isnan(v[p])])),
                         "missing": len(miss[k]),
                         "missing_positions": miss[k],
                         "impute_value_angstrom": round(float(np.nanmax(list(v.values()))), 2)}
                     for k, v in variants.items()},
        "pairwise_agreement": pairs_out,
        "ligand_reference_set_audit": lig_audit,
    }
    log_path = os.path.join(ROOT, "logs", "option_ab_distance_agreement.json")
    with open(log_path, "w") as f:
        json.dump(out, f, indent=1)
    print("")
    print("Wrote %s" % os.path.relpath(log_path, ROOT))

    print("")
    print("=" * 78)
    print("KEY CONTRASTS")
    print("  A_as_executed vs B_shell19  -> pure COORDINATE provenance (shell held fixed)")
    print("  B_as_executed vs B_shell19  -> pure REFERENCE-SET effect (coords held fixed)")
    print("  B_shell19     vs C_shell19  -> experimental coordinate replicate noise")
    print("=" * 78)
    print("")
    print("NO LABELS, FOLD-CHANGES OR PHENOTYPES WERE READ. STOP.")


if __name__ == "__main__":
    main()
