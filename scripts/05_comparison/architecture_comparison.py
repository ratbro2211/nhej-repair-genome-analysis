import pandas as pd

# ==========================
# LOAD FILES
# ==========================
hmm = pd.read_csv("hmm_subset_925.csv")
published = pd.read_csv("published_subset_925.csv")
interpro = pd.read_csv("interpro_architecture_subset_925.csv")

print("Files loaded")

# ==========================
# FIX COLUMN NAMES
# ==========================
hmm = hmm.rename(columns={"Architecture": "HMM"})
published = published.rename(columns={"Architecture": "Published"})
interpro = interpro.rename(columns={"Architecture": "InterPro"})

# ==========================
# FIX GENOME FORMAT (INTERPRO)
# ==========================
def fix_interpro(x):
    if pd.isna(x):
        return None
    if "GCF_" in x:
        return "GCF_" + x.split("GCF_")[1].split(".")[0]
    return None

interpro["genome_clean"] = interpro["tree_leaf"].apply(fix_interpro)

# ==========================
# KEEP REQUIRED COLUMNS
# ==========================
hmm = hmm[["genome_clean", "HMM"]]
published = published[["genome_clean", "Published"]]
interpro = interpro[["genome_clean", "InterPro"]]

# ==========================
# MERGE ALL THREE
# ==========================
df = hmm.merge(published, on="genome_clean")
df = df.merge(interpro, on="genome_clean")

print("\nMerged genomes:", len(df))

# ==========================
# DEFINE ALL 11 ARCHITECTURES
# ==========================
architectures = [
    "No-LigD-domains",
    "Ku-only",
    "LigD-only-all-same-prot",
    "LigD-only-all-diff-prot",
    "LigD-only-lig-no-pol",
    "LigD-only-lig-pol-diff-prot",
    "Yes-NHEJ-all-same-prot",
    "Yes-NHEJ-all-diff-prot",
    "Yes-NHEJ-lig-pol-diff-prot",
    "Yes-NHEJ-lig-pol-same_prot",
    "Yes-NHEJ-lig-no-pol"
]

# ==========================
# COUNT
# ==========================
hmm_counts = df["HMM"].value_counts()
pub_counts = df["Published"].value_counts()
int_counts = df["InterPro"].value_counts()

# ==========================
# BUILD FINAL TABLE
# ==========================
results = []

for arch in architectures:
    h = hmm_counts.get(arch, 0)
    p = pub_counts.get(arch, 0)
    i = int_counts.get(arch, 0)

    results.append([
        arch,
        p,
        h,
        i,
        h - p,
        i - p
    ])

comparison_df = pd.DataFrame(results, columns=[
    "Architecture",
    "Published",
    "HMM",
    "InterPro",
    "Diff (HMM-Published)",
    "Diff (InterPro-Published)"
])

# ==========================
# SAVE OUTPUTS
# ==========================
comparison_df.to_csv("Final_3way_comparison_925.csv", index=False)
df.to_csv("Genome_level_3way_comparison_925.csv", index=False)

# ==========================
# PRINT OUTPUT
# ==========================
print("\nFinal 3-way Comparison:\n")
print(comparison_df)

print("\nSaved files:")
print(" - Final_3way_comparison_925.csv")
print(" - Genome_level_3way_comparison_925.csv")
