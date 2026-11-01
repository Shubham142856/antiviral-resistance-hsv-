import pathlib
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from download_and_audit_structures import get_hsv1_reference_sequences, parse_chain_ca

ref_seqs = get_hsv1_reference_sequences()
parser = MMCIFParser(QUIET=True)

# 8V1Q
struct_8v1q = parser.get_structure('8V1Q', str(pathlib.Path(__file__).resolve().parent.parent / 'data/structures/8v1q.cif'))
chain_a_8v1q = parse_chain_ca(struct_8v1q, 'A')
ref_ul30 = ref_seqs['UL30']

print("--- 8V1Q Chain A Mismatches with NC_001806.2 UL30 ---")
for r in chain_a_8v1q:
    pos = r['resseq']
    if 1 <= pos <= len(ref_ul30):
        ref_aa = ref_ul30[pos - 1]
        if ref_aa != r['aa']:
            print(f"Pos {pos}: Structure={r['aa']} vs Reference={ref_aa}")

# 1OF1
struct_1of1 = parser.get_structure('1OF1', str(pathlib.Path(__file__).resolve().parent.parent / 'data/structures/1of1.cif'))
chain_b_1of1 = parse_chain_ca(struct_1of1, 'B')
ref_ul23 = ref_seqs['UL23']

print("\n--- 1OF1 Chain B Mismatches with NC_001806.2 UL23 ---")
for r in chain_b_1of1:
    pos = r['resseq']
    if 1 <= pos <= len(ref_ul23):
        ref_aa = ref_ul23[pos - 1]
        if ref_aa != r['aa']:
            print(f"Pos {pos}: Structure={r['aa']} vs Reference={ref_aa}")

