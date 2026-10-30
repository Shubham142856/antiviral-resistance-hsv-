import sqlite3
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq

ref_dir = 'herpesdrg/inst/ref'
for db in sorted(glob.glob(f'{ref_dir}/*.sqlite')):
    name = os.path.basename(db)
    conn = sqlite3.connect(db)
    print(f'=== DB: {name} ===')
    tables = [r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()]
    print('Tables:', tables)
    for t in tables:
        rows = conn.execute(f"SELECT * FROM {t} LIMIT 5").fetchall()
        print(f' Table {t}: {len(rows)} sample rows')
        for r in rows[:2]:
            print('  ', r)

