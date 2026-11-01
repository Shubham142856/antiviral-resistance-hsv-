import pathlib
import numpy as np
from Bio.PDB import MMCIFParser
from download_and_audit_structures import parse_chain_ca

parser = MMCIFParser(QUIET=True)

s_8v1q = parser.get_structure('8V1Q', str(pathlib.Path(__file__).resolve().parent.parent / 'data/structures/8v1q.cif'))
s_8exx = parser.get_structure('8EXX', str(pathlib.Path(__file__).resolve().parent.parent / 'data/structures/8exx.cif'))
s_8v1t = parser.get_structure('8V1T', str(pathlib.Path(__file__).resolve().parent.parent / 'data/structures/8v1t.cif'))

ca_8v1q = {r['resseq']: np.array([r['x'], r['y'], r['z']]) for r in parse_chain_ca(s_8v1q, 'A')}
ca_8exx = {r['resseq']: np.array([r['x'], r['y'], r['z']]) for r in parse_chain_ca(s_8exx, 'A')}
ca_8v1t = {r['resseq']: np.array([r['x'], r['y'], r['z']]) for r in parse_chain_ca(s_8v1t, 'A')}

common_q_exx = set(ca_8v1q.keys()).intersection(set(ca_8exx.keys()))
dists_q_exx = [np.linalg.norm(ca_8v1q[k] - ca_8exx[k]) for k in common_q_exx]
rmsd_q_exx = np.sqrt(np.mean(np.square(dists_q_exx)))

common_q_v1t = set(ca_8v1q.keys()).intersection(set(ca_8v1t.keys()))
dists_q_v1t = [np.linalg.norm(ca_8v1q[k] - ca_8v1t[k]) for k in common_q_v1t]
rmsd_q_v1t = np.sqrt(np.mean(np.square(dists_q_v1t)))

print(f"8V1Q vs 8EXX: common residues = {len(common_q_exx)}, C-alpha RMSD = {rmsd_q_exx:.3f} A, median distance = {np.median(dists_q_exx):.3f} A")
print(f"8V1Q vs 8V1T: common residues = {len(common_q_v1t)}, C-alpha RMSD = {rmsd_q_v1t:.3f} A, median distance = {np.median(dists_q_v1t):.3f} A")

