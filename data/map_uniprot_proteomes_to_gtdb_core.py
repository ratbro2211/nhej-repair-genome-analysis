import pandas as pd
import re

UNIPROT_TSV = "proteomes_taxonomy_id_2_AND_proteome_ty_2025_12_26.tsv"
GTDB_TSV    = "bac120_metadata_r226.tsv"
OUTPUT_TSV  = "uniprot_proteomes_with_GTDB_R226.tsv"

def core_accession(acc):
    if pd.isna(acc):
        return None
    acc = acc.strip()
    # Remove RS_/GB_ prefix + GCA_/GCF_
    return re.sub(r"^(RS_|GB_)?GC[AF]_", "", acc)

print("🔄 Loading UniProt proteome table...")
uni = pd.read_csv(UNIPROT_TSV, sep="\t", dtype=str)

print("🔍 UniProt columns:")
print(list(uni.columns))

# ---- IMPORTANT ----
# Adjust this column name if slightly different
ASSEMBLY_COL = "Genome assembly ID"

if ASSEMBLY_COL not in uni.columns:
    raise ValueError(f"Column '{ASSEMBLY_COL}' not found in UniProt file")

uni["assembly_core"] = uni[ASSEMBLY_COL].apply(core_accession)

print("🔄 Loading GTDB metadata...")
gtdb = pd.read_csv(GTDB_TSV, sep="\t", dtype=str)
gtdb["assembly_core"] = gtdb["accession"].apply(core_accession)

# Keep only relevant GTDB columns
gtdb_keep = gtdb[
    [
        "assembly_core",
        "accession",
        "gtdb_taxonomy",
        "gtdb_representative",
        "gtdb_genome_representative",
        "ncbi_organism_name",
        "ncbi_taxonomy",
        "ncbi_taxid",
        "ncbi_refseq_category",
        "ncbi_type_material_designation",
    ]
]

print("🔗 Mapping UniProt proteomes → GTDB (core accession match)...")
merged = uni.merge(
    gtdb_keep,
    on="assembly_core",
    how="left"
)

merged = merged.drop(columns=["assembly_core"])

merged.to_csv(OUTPUT_TSV, sep="\t", index=False)

print("✅ UniProt → GTDB mapping complete")
print(f"📁 Output file: {OUTPUT_TSV}")

