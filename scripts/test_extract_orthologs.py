import sqlite3
from Bio import SeqIO
from Bio.Seq import Seq

def get_gene_seq(ref, gene_name):
    db_path = f'herpesdrg/inst/ref/{ref}.sqlite'
    fasta_path = f'herpesdrg/inst/ref/{ref}.fasta'
    
    genome = str(SeqIO.read(fasta_path, 'fasta').seq)
    conn = sqlite3.connect(db_path)
    
    # Query transcript
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
        # try gene table
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
    
    print(f'Lookup {ref} - {gene_name}: {len(rows)} CDS segments found')
    if not rows:
        return None
    
    strand = rows[0][2]
    # assemble CDS
    cds_dna = ""
    for r in rows:
        c_start, c_end = r[5], r[6]
        # 1-based inclusive coordinates in TxDb
        seg = genome[c_start-1:c_end]
        cds_dna += seg
        
    s = Seq(cds_dna)
    if strand == '-':
        s = s.reverse_complement()
    prot = str(s.translate(to_stop=True))
    print(f'  Length DNA: {len(cds_dna)} bp, Protein: {len(prot)} aa. First 20: {prot[:20]}')
    return prot

# Check HSV-1
print("--- HSV-1 ---")
get_gene_seq('NC_001806.2', 'UL23')
get_gene_seq('NC_001806.2', 'UL30')

# Check HSV-2
print("--- HSV-2 ---")
get_gene_seq('NC_001798.2', 'UL23')
get_gene_seq('NC_001798.2', 'UL30')

# Check VZV (ORF36 = TK, ORF28 = DNA Pol)
print("--- VZV ---")
get_gene_seq('NC_001348.1', 'ORF36')
get_gene_seq('NC_001348.1', 'ORF28')

# Check HCMV (UL54 = DNA Pol)
print("--- HCMV ---")
get_gene_seq('NC_006273.2', 'UL54')

# Check HHV-6B (U38 = DNA Pol)
print("--- HHV-6B ---")
get_gene_seq('NC_000898.1', 'U38')

