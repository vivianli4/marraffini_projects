#!/usr/bin/env python3
import pandas as pd

# -------------------- CONFIG --------------------
INPUT_FILE = "ca4_foldseek_pdb100"
NODES_FILE = "nodes.dmp"
NAMES_FILE = "names.dmp"
OUTPUT_FILE = "filtered_results_ca4_prok.csv"

LDDT_CUTOFF = 0.2
# ------------------------------------------------


# -------------------- LOAD FOLDSEEK --------------------
def load_results(path):
    df = pd.read_csv(path, sep="\t", header=None)

    df = df.iloc[:, :9]

    df.columns = [
        "query", "target", "fident", "alnlen",
        "qcov", "tcov", "evalue", "bits", "lddt"
    ]

    num_cols = ["fident", "alnlen", "qcov", "tcov", "evalue", "bits", "lddt"]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


# -------------------- PARSE TARGET --------------------
def parse_target(df):
    df["pdb_id"] = df["target"].str.extract(r"^([0-9a-zA-Z]{4})")[0].str.lower()
    df["chain"] = df["target"].str.split("_").str[-1].str.upper()
    return df


# -------------------- LOAD NCBI TAXONOMY TREE --------------------
def load_ncbi_taxonomy(nodes_file):
    tax_ids = []
    parent_ids = []

    with open(nodes_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t|\t")

            try:
                tax_id = int(parts[0])
                parent_id = int(parts[1])

                tax_ids.append(tax_id)
                parent_ids.append(parent_id)

            except (ValueError, IndexError):
                continue  # skip malformed lines safely

    parent_map = dict(zip(tax_ids, parent_ids))
    return parent_map

# -------------------- LINEAGE RESOLVER --------------------
def get_lineage_fn(parent_map):
    def get_lineage(tax_id):
        lineage = set()
        current = tax_id

        while current in parent_map and current != parent_map[current]:
            lineage.add(current)
            current = parent_map[current]

        lineage.add(current)
        return lineage

    return get_lineage


# -------------------- CLASSIFY TAXONOMY --------------------
def annotate_superkingdom(df, get_lineage):
    BACTERIA = 2
    ARCHAEA = 2157
    VIRUSES = 10239
    EUKARYOTA = 2759

    def classify(tax_id):
        if pd.isna(tax_id):
            return "unknown"

        lineage = get_lineage(int(tax_id))

        if EUKARYOTA in lineage:
            return "eukaryote"
        if BACTERIA in lineage:
            return "bacteria"
        if ARCHAEA in lineage:
            return "archaea"
        if VIRUSES in lineage:
            return "virus"

        return "other"

    df["superkingdom"] = df["TAX_ID"].apply(classify)
    return df


# -------------------- MERGE TAXONOMY --------------------
def merge_taxonomy(df, sifts):
    df = df.merge(
        sifts,
        left_on=["pdb_id", "chain"],
        right_on=["PDB", "CHAIN"],
        how="left"
    )

    df["TAX_ID"] = pd.to_numeric(df["TAX_ID"], errors="coerce")
    return df


# -------------------- FINAL FILTER (CORRECT VERSION) --------------------
def filter_final(df):
    name = df["SCIENTIFIC_NAME"].str.lower().fillna("")

    is_phage = name.str.contains("phage|bacteriophage")

    return df[
        (df["superkingdom"].isin(["bacteria", "archaea"])) |
        (is_phage)
    ]


# -------------------- QUALITY FILTER --------------------
def apply_quality_filters(df):
    mask = pd.Series(True, index=df.index)

    if LDDT_CUTOFF is not None:
        mask &= df["lddt"] > LDDT_CUTOFF

    return df[mask]


# -------------------- FINAL OUTPUT --------------------
def finalize(df):
    df = df.sort_values(["query", "lddt"], ascending=[True, False])
    df = df.drop_duplicates(["query", "target"])

    return df[[
        "query",
        "target",
        "evalue",
        "lddt",
        "pdb_id",
        "chain",
        "SCIENTIFIC_NAME"
    ]]


# -------------------- MAIN --------------------
def main():
    print("Loading Foldseek...")
    df = load_results(INPUT_FILE)

    print("Parsing targets...")
    df = parse_target(df)

    print("Loading taxonomy tree...")
    parent_map = load_ncbi_taxonomy(NODES_FILE)
    get_lineage = get_lineage_fn(parent_map)

    print("Loading SIFTS...")
    sifts = pd.read_csv("pdb_chain_taxonomy.tsv", sep="\t", dtype=str, comment="#")
    sifts.columns = [c.strip().upper() for c in sifts.columns]

    print("Merging taxonomy...")
    df = merge_taxonomy(df, sifts)

    print("Annotating lineage...")
    df = annotate_superkingdom(df, get_lineage)

    print("Filtering bacteria + phage...")
    df = filter_final(df)

    print("Applying quality filters...")
    df = apply_quality_filters(df)

    print("Finalizing...")
    df = finalize(df)

    print(f"Final rows: {len(df)}")

    df.to_csv(OUTPUT_FILE, index=False)
    print(f"Saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()