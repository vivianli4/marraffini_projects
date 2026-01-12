#!/usr/bin/env python3

import pandas as pd
import subprocess
import os

# ----- Step 1: Read CSV -----
csv_file = "all_designs_metrics_ca4.csv"  # replace with your CSV filename
df = pd.read_csv(csv_file)

# ----- Step 2: Convert CSV to FASTA -----
fasta_file = "ca4_binders.fasta"
with open(fasta_file, "w") as f:
    for index, row in df.iterrows():
        f.write(f">{row['id']}\n{row['designed_sequence']}\n")

print(f"FASTA file created: {fasta_file}")

# ----- Step 3: Run BLASTP -----
# Make sure you have a BLAST database, e.g., 'swissprot'
blast_db = "swissprot"   # change to your BLAST database name
output_file = "blast_results_ca4.txt"

# Run BLASTP
blast_command = [
    "blastp",
    "-query", fasta_file,
    "-db", blast_db,
    "-out", output_file,
    "-outfmt", "6"  # tabular output
]

subprocess.run(blast_command)
print(f"BLASTP results saved to: {output_file}")