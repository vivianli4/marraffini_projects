#!/usr/bin/env python3
import pandas as pd

#change the name of the file read 

df = pd.read_csv(
    "swissprot_hits.tsv",
    sep="\t",
    header=None,
    names=[
        "query",
        "target",
        "evalue",
        "bits",
        "qcov",
        "tcov",
        "alnlen",
    ],
    dtype=str
)

# Convert numeric columns AFTER loading
df["evalue"] = df["evalue"].astype(float)
df["bits"] = df["bits"].astype(float)


# Filter by E-value
df = df[df["evalue"] < 1e-4]

# Extract UniProt ID
df["uniprot"] = df["target"].str.extract(r"AF-(.*?)-F1")

df = df.dropna(subset=["uniprot"])

df[["uniprot"]].drop_duplicates().to_csv(
    "uniprot_ids_ca4.txt",
    index=False,
    header=False
)

print(f"Saved {df['uniprot'].nunique()} UniProt IDs")
