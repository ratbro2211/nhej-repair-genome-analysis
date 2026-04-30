#!/usr/bin/env python3

import pandas as pd
from Bio import Phylo

TREE_FILE = "bac120_r226.tree"
TSV_FILE = "uniprot_tree_placement.tsv"
OUT_TREE = "bac120_r226_uniprot_pruned.tree"

print("🔄 Loading GTDB tree...")
tree = Phylo.read(TREE_FILE, "newick")
tree_leaves = {t.name for t in tree.get_terminals() if t.name}
print(f"  Tree leaves: {len(tree_leaves):,}")

print("\n🔄 Loading UniProt placement TSV...")
df = pd.read_csv(TSV_FILE, sep="\t")

required = {"accession", "gtdb_genome_representative", "placement_status"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing columns: {missing}")

# ---------------------------------------------------
# RESOLVE ALL ROWS TO TREE LEAVES (CORE FIX)
# ---------------------------------------------------
resolved_leaves = set()

for _, row in df.iterrows():
    status = row["placement_status"]

    if status == "direct_match":
        resolved_leaves.add(row["accession"])

    elif status == "mapped_to_representative":
        resolved_leaves.add(row["gtdb_genome_representative"])

    # unplaceable → ignored

print("\n📊 Placement summary")
print(f"  Direct match leaves           : {(df['placement_status']=='direct_match').sum():,}")
print(f"  Representative-mapped leaves  : {(df['placement_status']=='mapped_to_representative').sum():,}")
print(f"  Unique resolved tree leaves   : {len(resolved_leaves):,}")
print(f"  TSV rows                      : {len(df):,}")

# ---------------------------------------------------
# PRUNE TREE
# ---------------------------------------------------
to_remove = tree_leaves - resolved_leaves
print(f"\n✂️  Pruning {len(to_remove):,} unused GTDB leaves...")

for leaf in to_remove:
    try:
        tree.prune(leaf)
    except ValueError:
        pass

Phylo.write(tree, OUT_TREE, "newick")

print("\n✅ DONE")
print(f"🌲 Final tree leaves : {len(tree.get_terminals()):,}")
print(f"📁 Saved as          : {OUT_TREE}")

