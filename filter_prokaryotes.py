#!/usr/bin/env python3

#April 8, 2026 written with chatbot

import pandas as pd
import requests

INPUT_FILE = "master_table_ca6.tsv"
OUTPUT_FILE = "filtered_ca6_bacterial.tsv"

# Load table
df = pd.read_csv(INPUT_FILE, sep="\t")

# Get unique UniProt IDs (avoid redundant API calls)
unique_ids = df["uniprot"].dropna().unique()

print(f"Total unique UniProt IDs: {len(unique_ids)}")

# Function to query UniProt
def get_taxonomy(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return None, None
        
        data = r.json()
        lineage = data["organism"]["lineage"]
        
        # Extract superkingdom (first entry)
        superkingdom = lineage[0] if lineage else None
        
        return lineage, superkingdom
    
    except Exception as e:
        return None, None

# Build mapping: UniProt ID → taxonomy
taxonomy_dict = {}

for i, uid in enumerate(unique_ids):
    if i % 50 == 0:
        print(f"Processing {i}/{len(unique_ids)}")
    
    lineage, superkingdom = get_taxonomy(uid)
    
    taxonomy_dict[uid] = {
        "lineage": ";".join(lineage) if lineage else None,
        "superkingdom": superkingdom
    }
    

# Map taxonomy back to dataframe
df["lineage"] = df["uniprot"].map(lambda x: taxonomy_dict.get(x, {}).get("lineage"))
df["superkingdom"] = df["uniprot"].map(lambda x: taxonomy_dict.get(x, {}).get("superkingdom"))

# Filter: keep only Bacteria + Archaea
filtered_df = df[df["superkingdom"].isin(["Bacteria", "Archaea"])]

print(f"Original rows: {len(df)}")
print(f"Filtered rows: {len(filtered_df)}")

# Save
filtered_df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print(f"Filtered table saved to: {OUTPUT_FILE}")