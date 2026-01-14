#!/usr/bin/env python3
import pandas as pd
from Bio import Entrez
from pathlib import Path
import subprocess
import json
import time
import sys
import ssl
import certifi
import urllib.request

ssl._create_default_https_context = ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())


FOLDSEEK_FILE = "swissprot_hits.tsv"
EVALUE_CUTOFF = 1e-4

Entrez.email = "yli03@rockefeller.edu"
Entrez.api_key = None  # optional but recommended

SLEEP = 0.3  # NCBI-safe rate

df = pd.read_csv(
    FOLDSEEK_FILE,
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

print("Mapping UniProt RefSeq proteins...")

uniprots= Path('uniprot_ids_ca4.txt').read_text().splitlines()
refseq_proteins = set()

for up in uniprots:
    handle = Entrez.esearch(
        db="protein",
        term=f"{up}[Accession] AND refseq[filter]"
    )
    record = Entrez.read(handle)
    handle.close()

    if record["IdList"]:
        refseq_proteins.add(record["IdList"][0])

    time.sleep(SLEEP)

Path("refseq_proteins_ca4.txt").write_text("\n".join(sorted(refseq_proteins)))

print(f"Saved {len(refseq_proteins)} RefSeq protein accessions")

############################
# STEP 3: RefSeq → Assemblies (Datasets)
############################

print("Mapping RefSeq proteins → genome assemblies (Datasets)...")

if not refseq_proteins:
    print("ERROR: No RefSeq proteins found. Stopping.")
    sys.exit(1)

cmd = [
    "datasets", "summary", "protein", "accession",
    "--inputfile", "refseq_proteins.txt",
    "--report", "genome",
    "--as-json"
]

result = subprocess.run(
    cmd,
    capture_output=True,
    text=True,
    check=True
)

data = json.loads(result.stdout)

assemblies = set()

for report in data.get("reports", []):
    asm = report.get("assembly_accession")
    if asm:
        assemblies.add(asm)

Path("assemblies_ca4.txt").write_text("\n".join(sorted(assemblies)))

print(f"Saved {len(assemblies)} genome assemblies")

print("Pipeline complete ✅")
