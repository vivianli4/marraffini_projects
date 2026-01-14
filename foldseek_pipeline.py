#!/usr/bin/env python3
import subprocess
import pandas as pd
from pathlib import Path

QUERIES_DIR = "final_100_designs_cam1"
DB = "af_swissprot"
OUT = "cam1_foldseek"
TMP = "tmp"
THREADS = 16


def run_foldseek():
    cmd = [
        "foldseek", "easy-search",
        QUERIES_DIR,
        DB,
        OUT,
        TMP,
        "--threads", str(THREADS),
        "--max-seqs", "1000",
        "--format-output",
        "query,target,fident,alnlen,qcov,tcov,evalue,bits",
    ]
    subprocess.run(cmd, check=True)


def load_results():
    cols = [
        "query",
        "target",
        "fident",
        "alnlen",
        "qcov",
        "tcov",
        "evalue",
        "bits",
    ]
    return pd.read_csv(OUT, sep="\t", names=cols)


def postprocess(df):
    # Extract UniProt from AFDB ID
    df["uniprot"] = df["target"].str.extract(r"AF-(.*?)-F1")

    # Example filters
    df = df[
        (df["evalue"] < 1e-5) &
        (df["qcov"] > 0.5)
    ]

    return df.sort_values(
        ["query", "bits"],
        ascending=[True, False]
    )


if __name__ == "__main__":
    run_foldseek()
    df = load_results()
    df = postprocess(df)

    df.to_csv("filtered_results.csv", index=False)
    print(df.head())
