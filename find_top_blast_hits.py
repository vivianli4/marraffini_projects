#!/usr/bin/env python3
import pandas as pd

# ---- Step 1: Load BLAST results ----
blast_file = "blast_results_cam1.txt"  # change this if your file has a different name

# Columns for BLAST tabular output
cols = ["query_id","subject_id","perc_identity","align_len","mismatches","gap_opens",
        "q_start","q_end","s_start","s_end","evalue","bit_score"]

# Load into pandas DataFrame
df = pd.read_csv(blast_file, sep="\t", names=cols)

# ---- Step 2: Filter by E-value < 0.1 ----
filtered = df[df.evalue < 0.1]

# ---- Step 3: Sort by query and E-value ----
filtered = filtered.sort_values(["query_id","evalue"])

filtered.to_csv("blast_results_filtered.txt", sep="\t", index=False)