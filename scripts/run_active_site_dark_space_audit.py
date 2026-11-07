"""
run_active_site_dark_space_audit.py
===================================
Rigorous quantification of the "Clinical Resistance Dark Space" across the
catalytic active-site shells of HSV-1 Thymidine Kinase (UL23) and DNA
Polymerase (UL30).

Quantifies:
  1. The ascertainment coverage of 3D catalytic contact shells (<= 5.0 A from ligand/cofactors).
  2. The biological mechanism explaining the coverage disparity:
     - UL23 (dispensable enzyme): 81.2% characterized (13/16), 18.8% dark space.
     - UL30 (essential enzyme): 26.9% characterized (7/26), 73.1% dark space.
  3. Prospective Clinical Surveillance Watchlist for the 19 uncharacterized UL30 residues
     and 3 uncharacterized UL23 residues, integrating cryo-EM ligand distances,
     multi-herpesvirus MSA conservation, and phenotypic risk tiering.

Outputs:
  - logs/active_site_dark_space_surveillance.json
"""

import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import json
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.Align import PairwiseAligner, substitution_matrices

def compute_msa_maps():
    # Define ortholog sequences directly from GenBank references
    from scripts.compute_msa_conservation_baseline import orthologs, align_to_ref
    
    entropy_map = {}
    freq_conserved_map = {}

    for gene, o_dict in orthologs.items():
        ref_seq = o_dict["HSV-1"]
        L = len(ref_seq)
        pos_residues = {p: [ref_seq[p-1]] for p in range(1, L + 1)}

        for o_name, q_seq in o_dict.items():
            if o_name == "HSV-1":
                continue
            mapping = align_to_ref(ref_seq, q_seq)
            for p in range(1, L + 1):
                ch = mapping.get(p, "-")
                pos_residues[p].append(ch)

        for p in range(1, L + 1):
            res_list = pos_residues[p]
            valid_aas = [ch for ch in res_list if ch != "-"]
            wt_aa = ref_seq[p-1]

            if valid_aas:
                f_cons = sum(1 for ch in valid_aas if ch == wt_aa) / len(valid_aas)
            else:
                f_cons = 0.0
            freq_conserved_map[(gene, p)] = f_cons

            if len(valid_aas) <= 1:
                h = 0.0
            else:
                counts = pd.Series(valid_aas).value_counts()
                probs = counts / counts.sum()
                h = -float(np.sum(probs * np.log2(probs)))
            entropy_map[(gene, p)] = h

    return orthologs, freq_conserved_map, entropy_map


def main():
    print("=" * 80)
    print("RUNNING ACTIVE-SITE DARK-SPACE SURVEILLANCE MAPPING")
    print("=" * 80)

    # 1. Load clinical cohort positions
    feat_path = os.path.join(ROOT, "data", "processed", "where_what_features.tsv")
    df_hsv = pd.read_csv(feat_path, sep="\t")
    ul23_clin = set(df_hsv[df_hsv["gene"] == "UL23"]["pos"].dropna().astype(int).unique())
    ul30_clin = set(df_hsv[df_hsv["gene"] == "UL30"]["pos"].dropna().astype(int).unique())

    # 2. Compute MSA conservation and entropy profiles
    orthologs, conserv_data, entropy_data = compute_msa_maps()
    ref_ul30 = orthologs["UL30"]["HSV-1"]
    ref_ul23 = orthologs["UL23"]["HSV-1"]

    cif_parser = MMCIFParser(QUIET=True)

    # ----------------------------------------------------
    # UL30 Active Site Shell (8EXX foscarnet, 8V1T aciclovir)
    # ----------------------------------------------------
    print("Parsing cryo-EM UL30 structures (8EXX foscarnet, 8V1T aciclovir)...")
    s_8exx = cif_parser.get_structure("8EXX", os.path.join(ROOT, "data", "structures", "8exx.cif"))
    m_8exx = list(s_8exx.get_models())[0]
    ppf_atoms = []
    for res in m_8exx.get_residues():
        if res.get_resname() in ["PPF", "MG"]:
            for a in res.get_atoms():
                ppf_atoms.append(a.get_coord())
    ppf_coords = np.array(ppf_atoms)

    s_8v1t = cif_parser.get_structure("8V1T", os.path.join(ROOT, "data", "structures", "8v1t.cif"))
    m_8v1t = list(s_8v1t.get_models())[0]
    avp_atoms = []
    for res in m_8v1t.get_residues():
        if res.get_resname() in ["AVP", "ACV"]:
            for a in res.get_atoms():
                avp_atoms.append(a.get_coord())
    avp_coords = np.array(avp_atoms)

    # Map all residues in chain A of 8EXX / 8V1T within 5.0 A of ligand
    ul30_as_all = {}
    for res in m_8exx["A"].get_residues():
        if res.id[0] == " ":
            pos = res.id[1]
            coords = np.array([a.get_coord() for a in res.get_atoms()])
            d_ppf = float(np.min(np.linalg.norm(ppf_coords - coords[:, None, :], axis=2)))
            ul30_as_all[pos] = {"dist_ppf": d_ppf, "dist_avp": 999.0}

    for res in m_8v1t["A"].get_residues():
        if res.id[0] == " ":
            pos = res.id[1]
            coords = np.array([a.get_coord() for a in res.get_atoms()])
            d_avp = float(np.min(np.linalg.norm(avp_coords - coords[:, None, :], axis=2)))
            if pos in ul30_as_all:
                ul30_as_all[pos]["dist_avp"] = d_avp
            else:
                ul30_as_all[pos] = {"dist_ppf": 999.0, "dist_avp": d_avp}

    # Filter to <= 5.0 A
    ul30_active_positions = {}
    for p, d in ul30_as_all.items():
        min_d = min(d["dist_ppf"], d["dist_avp"])
        if min_d <= 5.0:
            ul30_active_positions[p] = {
                "min_dist": min_d,
                "dist_ppf": d["dist_ppf"],
                "dist_avp": d["dist_avp"]
            }

    # ----------------------------------------------------
    # UL23 Active Site Shell (1OF1 ganciclovir/thymidine)
    # ----------------------------------------------------
    print("Parsing crystal structure UL23 (1OF1)...")
    s_1of1 = cif_parser.get_structure("1OF1", os.path.join(ROOT, "data", "structures", "1of1.cif"))
    m_1of1 = list(s_1of1.get_models())[0]
    sct_atoms = []
    for res in m_1of1.get_residues():
        if res.get_resname() in ["SCT", "THM"]:
            for a in res.get_atoms():
                sct_atoms.append(a.get_coord())
    sct_coords = np.array(sct_atoms)

    ul23_active_positions = {}
    for res in m_1of1["A"].get_residues():
        if res.id[0] == " ":
            pos = res.id[1]
            coords = np.array([a.get_coord() for a in res.get_atoms()])
            d_sct = float(np.min(np.linalg.norm(sct_coords - coords[:, None, :], axis=2)))
            if d_sct <= 5.0:
                ul23_active_positions[pos] = {
                    "min_dist": d_sct,
                    "dist_ligand": d_sct
                }

    # Build Surveillance Watchlist
    ul30_watchlist = []
    for p in sorted(ul30_active_positions.keys()):
        is_clin = (p in ul30_clin)
        wt = ref_ul30[p - 1]
        cons = conserv_data.get(("UL30", p), 0.0)
        ent = entropy_data.get(("UL30", p), 0.0)
        d_info = ul30_active_positions[p]

        # Biological risk stratification
        if p in [717, 888]:
            risk_tier = "LETHAL_CATALYTIC_CORE"
            mechanism = "Catalytic triad carboxylate; essential for coordination of Mg2+ and phosphoryl transfer. Mutations ablate replication."
        elif d_info["dist_ppf"] <= 3.5:
            risk_tier = "HIGH_FOSCARNET_POCKET"
            mechanism = "Direct pyrophosphate-binding shell; high barrier to non-synonymous mutation, selective for pyrophosphate analogs."
        elif d_info["dist_avp"] <= 3.5:
            risk_tier = "HIGH_NUCLEOTIDE_POCKET"
            mechanism = "Direct dNTP/aciclovir triphosphate-binding pocket; potential discriminative resistance if polymerase kinetics preserved."
        else:
            risk_tier = "MODERATE_CONTACT_SHELL"
            mechanism = "Second-shell catalytic contact; modulates active center conformation and fidelity."

        ul30_watchlist.append({
            "gene": "UL30",
            "pos": p,
            "wt": wt,
            "min_dist_angstrom": d_info["min_dist"],
            "dist_ppf_angstrom": d_info["dist_ppf"],
            "dist_avp_angstrom": d_info["dist_avp"],
            "msa_conservation": float(cons),
            "shannon_entropy": float(ent),
            "status": "CLINICAL_CHARACTERIZED" if is_clin else "DARK_SPACE_WATCHLIST",
            "risk_tier": risk_tier,
            "biological_mechanism": mechanism
        })

    ul23_watchlist = []
    for p in sorted(ul23_active_positions.keys()):
        is_clin = (p in ul23_clin)
        wt = ref_ul23[p - 1]
        cons = conserv_data.get(("UL23", p), 0.0)
        ent = entropy_data.get(("UL23", p), 0.0)
        d_info = ul23_active_positions[p]

        if p in [101, 132]:
            risk_tier = "HIGH_AROMATIC_STACKING"
            mechanism = "Aromatic sandwiching of purine ring; missense variants likely induce substantial loss of acyclovir phosphorylation."
        elif p == 97:
            risk_tier = "MODERATE_HYDROPHOBIC_FLOOR"
            mechanism = "Hydrophobic substrate binding floor; subtle steric accommodation."
        else:
            risk_tier = "CHARACTERIZED_RESISTANCE_HOTSPOT"
            mechanism = "Well-characterized clinical drug-resistance or kinase inactivation position."

        ul23_watchlist.append({
            "gene": "UL23",
            "pos": p,
            "wt": wt,
            "min_dist_angstrom": d_info["min_dist"],
            "dist_ligand_angstrom": d_info["dist_ligand"],
            "msa_conservation": float(cons),
            "shannon_entropy": float(ent),
            "status": "CLINICAL_CHARACTERIZED" if is_clin else "DARK_SPACE_WATCHLIST",
            "risk_tier": risk_tier,
            "biological_mechanism": mechanism
        })

    # Summary statistics
    ul30_total = len(ul30_watchlist)
    ul30_dark = len([r for r in ul30_watchlist if r["status"] == "DARK_SPACE_WATCHLIST"])
    ul30_char = ul30_total - ul30_dark

    ul23_total = len(ul23_watchlist)
    ul23_dark = len([r for r in ul23_watchlist if r["status"] == "DARK_SPACE_WATCHLIST"])
    ul23_char = ul23_total - ul23_dark

    output_data = {
        "summary": {
            "UL30_catalytic_shell_residues": ul30_total,
            "UL30_clinically_characterized": ul30_char,
            "UL30_characterized_pct": float(ul30_char / ul30_total * 100),
            "UL30_dark_space_residues": ul30_dark,
            "UL30_dark_space_pct": float(ul30_dark / ul30_total * 100),
            "UL23_catalytic_shell_residues": ul23_total,
            "UL23_clinically_characterized": ul23_char,
            "UL23_characterized_pct": float(ul23_char / ul23_total * 100),
            "UL23_dark_space_residues": ul23_dark,
            "UL23_dark_space_pct": float(ul23_dark / ul23_total * 100),
            "biological_insight": (
                "Dispensable viral kinase (UL23) permits viral replication upon complete enzyme "
                "inactivation, enabling 81.2% saturation of the active pocket by clinical resistance. "
                "In contrast, essential replicative polymerase (UL30) imposes stringent viability "
                "constraints, leaving 73.1% of catalytic contact residues in unobserved dark space."
            )
        },
        "ul30_catalytic_shell": ul30_watchlist,
        "ul23_catalytic_shell": ul23_watchlist
    }

    out_file = os.path.join(ROOT, "logs", "active_site_dark_space_surveillance.json")
    with open(out_file, "w") as f:
        json.dump(output_data, f, indent=2)
    print(f"Saved active-site dark space audit to: {out_file}")

    print("\n" + "=" * 95)
    print(f"| {'Target':<10} | {'Shell (<=5A)':<14} | {'Characterized':<16} | {'Dark Space':<16} | {'Dark Space Pct':<16} |")
    print("|" + "-" * 12 + "|" + "-" * 16 + "|" + "-" * 18 + "|" + "-" * 18 + "|" + "-" * 18 + "|")
    print(f"| {'UL23 (TK)':<10} | {ul23_total:<14} | {ul23_char:<16} | {ul23_dark:<16} | {ul23_dark/ul23_total*100:>14.1f}% |")
    print(f"| {'UL30 (Pol)':<10} | {ul30_total:<14} | {ul30_char:<16} | {ul30_dark:<16} | {ul30_dark/ul30_total*100:>14.1f}% |")
    print("=" * 95)

if __name__ == "__main__":
    main()
