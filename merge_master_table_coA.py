#!/usr/bin/env python3

import pandas as pd

# Input files
FOLDSEEK_FILE = "filtered_results_ca4.csv"
MAPPING_FILE = "full_mapping_ca4.tsv"
OUTPUT_FILE = "master_table_ca4.tsv"

# Load Foldseek results
foldseek = pd.read_csv(FOLDSEEK_FILE)
foldseek = foldseek[["query", "uniprot", "target", "evalue"]]

# Load UniProt → RefSeq → Genome mapping
mapping = pd.read_csv(MAPPING_FILE, sep="\t")
mapping = mapping.rename(columns={
    "UniProt": "uniprot",
    "Protein": "refseq",
    "Assembly": "genome_accession"
})
mapping = mapping[mapping["genome_accession"] != "NA"]  # drop missing genomes
mapping["refseq"] = mapping["refseq"].str.split(".").str[0]  # remove version

# Merge Foldseek hits with genome mapping
df = foldseek.merge(mapping, on="uniprot", how="left")
df = df.dropna(subset=["genome_accession"])  # remove rows that didn't map
df = df[["genome_accession", "query", "uniprot", "refseq", "evalue", "target"]]
df = df.sort_values(["genome_accession", "query", "evalue"])

# Save master table
df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print(f"Master table saved: {OUTPUT_FILE}")
print(f"Rows: {len(df)}, Unique genomes: {df['genome_accession'].nunique()}, Unique synthetic proteins: {df['query'].nunique()}")
