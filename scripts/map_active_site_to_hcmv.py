"""
map_active_site_to_hcmv.py
==========================
Maps experimental ligand-contact residues from HSV-1 UL30 cryo-EM structures:
  - 8EXX: holoenzyme bound to foscarnet (PPF)
  - 8V1T: holoenzyme bound to acyclovir triphosphate (AVP)
onto HCMV UL54 reference sequence (NC_006273.2, 1242 aa) via global BLOSUM62 alignment.
"""

import os
import sqlite3
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def get_gene_seq(ref, gene_name):
    db_path = os.path.join(ROOT, "herpesdrg", "inst", "ref", f"{ref}.sqlite")
    fasta_path = os.path.join(ROOT, "herpesdrg", "inst", "ref", f"{ref}.fasta")
    genome = str(SeqIO.read(fasta_path, 'fasta').seq)
    conn = sqlite3.connect(db_path)
    query = """
    SELECT t._tx_id, t.tx_name, t.tx_strand, t.tx_start, t.tx_end, c.cds_start, c.cds_end
    FROM transcript t
    JOIN splicing s ON t._tx_id = s._tx_id
    JOIN cds c ON s._cds_id = c._cds_id
    WHERE t.tx_name = ? OR t.tx_name LIKE ?
    ORDER BY s.exon_rank
    """
    rows = conn.execute(query, (gene_name, f'%{gene_name}%')).fetchall()
    if not rows:
        query2 = """
        SELECT t._tx_id, t.tx_name, t.tx_strand, t.tx_start, t.tx_end, c.cds_start, c.cds_end
        FROM gene g
        JOIN transcript t ON g._tx_id = t._tx_id
        JOIN splicing s ON t._tx_id = s._tx_id
        JOIN cds c ON s._cds_id = c._cds_id
        WHERE g.gene_id = ? OR g.gene_id LIKE ?
        ORDER BY s.exon_rank
        """
        rows = conn.execute(query2, (gene_name, f'%{gene_name}%')).fetchall()
    if not rows:
        raise ValueError(f"Could not find CDS for {ref} {gene_name}")
    strand = rows[0][2]
    cds_dna = "".join([genome[r[5]-1:r[6]] for r in rows])
    s = Seq(cds_dna)
    if strand == '-':
        s = s.reverse_complement()
    prot = str(s.translate(to_stop=True))
    return prot

seq_hsv = get_gene_seq('NC_001806.2', 'UL30')
seq_hcmv = get_gene_seq('NC_006273.2', 'UL54')

print(f'HSV-1 UL30 length: {len(seq_hsv)} aa')
print(f'HCMV UL54 length: {len(seq_hcmv)} aa')

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10.0
aligner.extend_gap_score = -0.5
aligner.mode = 'global'

aln = aligner.align(seq_hsv, seq_hcmv)[0]
ref_aln = aln[0]
qry_aln = aln[1]

hsv_idx = 0
hcmv_idx = 0
hsv_to_hcmv = {}

for r_ch, q_ch in zip(ref_aln, qry_aln):
    if r_ch != '-':
        hsv_idx += 1
    if q_ch != '-':
        hcmv_idx += 1
    if r_ch != '-' and q_ch != '-':
        hsv_to_hcmv[hsv_idx] = (hcmv_idx, r_ch, q_ch)

# 19 contact residues experimentally identified in 8EXX (FOS) and 8V1T (ACV-TP) (cutoff <= 5.0 A)
contact_positions = [
    717, 718, 719, 720, 721, 722, 723, # Region II / Motif A
    785, 789,                          # Motif B
    811, 812, 815, 816, 818, 819,      # Motif B
    887, 888,                          # Motif C (catalytic triad)
    925, 927                           # Palm C-term
]

mapping_results = []
mapped_hcmv_positions = []

print('\nMapping of HSV-1 UL30 Experimental Drug Contacts to HCMV UL54:')
for p in contact_positions:
    if p in hsv_to_hcmv:
        hp, r_ch, q_ch = hsv_to_hcmv[p]
        mapped_hcmv_positions.append(hp)
        mapping_results.append({
            'hsv1_pos': p,
            'hsv1_wt': r_ch,
            'hcmv_pos': hp,
            'hcmv_wt': q_ch,
            'conserved': (r_ch == q_ch)
        })
        print(f"  HSV-1 {r_ch}{p:4d}  -->  HCMV UL54 {q_ch}{hp:4d}  {'[CONSERVED]' if r_ch == q_ch else '[SUBSTITUTED]'}")

print('\nSummary of Mapped HCMV UL54 Active-Site Residues (1-based indices):')
print(sorted(mapped_hcmv_positions))

with open('logs/hcmv_ul54_active_site_mapping.json', 'w') as f:
    json.dump({
        'source_experimental_pdbs': ['8EXX (FOS, PPF)', '8V1T (ACV-TP, AVP)'],
        'distance_cutoff_angstrom': 5.0,
        'mapped_residues': mapping_results,
        'hcmv_active_site_positions': sorted(mapped_hcmv_positions)
    }, f, indent=2)

print('\nLogged mapping to logs/hcmv_ul54_active_site_mapping.json')
