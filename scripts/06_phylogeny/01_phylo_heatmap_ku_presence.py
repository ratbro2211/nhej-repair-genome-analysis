"""
Phylogenetic heatmap aligned to GTDB bac120 R226 tree.

Builds 12 binary domain-presence columns (KU/LIG/POL/PE x PUB/HMM/IPR)
plus 2 categorical taxonomy tracks (Phylum, Class) extracted from the
GTDB taxonomy column.  All 14 tracks are drawn as discrete blocks
aligned to the pruned tree leaves.
"""

import os, re, sys, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_hex
from Bio import Phylo
from ete3 import Tree as EteTree

# -- paths -----------------------------------------------------------------
HMM_F   = f"{BASE}/hmm_subset_925.csv"
IPR_F   = f"{BASE}/interpro_architecture_final.csv"
PUB_F   = f"{BASE}/published_subset_925.csv"
GTDB_F  = f"{BASE}/uniprot_proteomes_with_GTDB_R226.tsv"
TREE_F  = f"{BASE}/bac120_r226.tree"
OUT_PDF = f"{BASE}/figures/phylo_heatmap.pdf"
OUT_PNG = f"{BASE}/figures/phylo_heatmap.png"
OUT_CSV = f"{BASE}/results/binary_matrix.csv"
OUT_TAX = f"{BASE}/results/taxonomy_table.csv"
os.makedirs(f"{BASE}/figures", exist_ok=True)
os.makedirs(f"{BASE}/results", exist_ok=True)

np.random.seed(42)

# -- helper: normalise genome ID -------------------------------------------
def norm_id(raw: str) -> str:
    """Remove RS_/GB_ prefix and trailing version (.N)."""
    s = str(raw).strip()
    s = re.sub(r'^(RS_|GB_)', '', s)
    s = re.sub(r'\.\d+$', '', s)
    return s


# -- 1. load datasets ------------------------------------------------------
print("=== 1. Loading datasets ===")

hmm  = pd.read_csv(HMM_F)
ipr  = pd.read_csv(IPR_F)
pub  = pd.read_csv(PUB_F)
gtdb = pd.read_csv(GTDB_F, sep="\t", usecols=["accession", "gtdb_taxonomy"])

# normalise IDs
hmm["gid"]  = hmm["genome_clean"].apply(norm_id)
ipr["gid"]  = ipr["tree_leaf"].apply(norm_id)
pub["gid"]  = pub["genome_clean"].apply(norm_id)
gtdb["gid"] = gtdb["accession"].apply(norm_id)

print(f"  HMM  : {len(hmm)} rows, unique IDs: {hmm['gid'].nunique()}")
print(f"  IPR  : {len(ipr)} rows, unique IDs: {ipr['gid'].nunique()}")
print(f"  PUB  : {len(pub)} rows, unique IDs: {pub['gid'].nunique()}")
print(f"  GTDB : {len(gtdb)} rows, unique IDs: {gtdb['gid'].nunique()}")

# keep only one row per genome (in case of duplicates)
hmm  = hmm.drop_duplicates("gid").set_index("gid")
ipr  = ipr.drop_duplicates("gid").set_index("gid")
pub  = pub.drop_duplicates("gid").set_index("gid")
gtdb = gtdb.drop_duplicates("gid").set_index("gid")


# -- 2. parse tree & get leaf IDs ------------------------------------------
print("\n=== 2. Reading tree ===")
t0 = time.time()
ete_tree = EteTree(TREE_F, format=1, quoted_node_names=True)
print(f"  Full tree parsed in {time.time()-t0:.1f}s")
tree_ids_raw = {leaf.name for leaf in ete_tree.get_leaves()}
tree_ids = {norm_id(n) for n in tree_ids_raw}
print(f"  Tree leaves (raw): {len(tree_ids_raw)}  (normalised unique): {len(tree_ids)}")

# build reverse map: normalised_id -> original_tree_name
norm2raw = {}
for name in tree_ids_raw:
    n = norm_id(name)
    norm2raw[n] = name


# -- 3. intersection -------------------------------------------------------
print("\n=== 3. Computing intersection ===")
sets = {
    "HMM":  set(hmm.index),
    "IPR":  set(ipr.index),
    "PUB":  set(pub.index),
    "GTDB": set(gtdb.index),
    "TREE": tree_ids,
}
for name, s in sets.items():
    print(f"  {name}: {len(s)}")

intersection = (set(hmm.index) & set(ipr.index) & set(pub.index)
                & set(gtdb.index) & set(tree_ids))
print(f"\n  Intersection (HMM & IPR & PUB & GTDB & TREE): {len(intersection)}")

if len(intersection) == 0:
    print("  WARNING: empty intersection - check ID normalisation")
    sys.exit(1)

genome_list = sorted(intersection)
print(f"  Using {len(genome_list)} genomes")


# -- 4. prune tree ---------------------------------------------------------
print("\n=== 4. Pruning tree ===")
keep_raw = [norm2raw[g] for g in genome_list if g in norm2raw]
print(f"  Pruning to {len(keep_raw)} leaves (ETE3)...")
t0 = time.time()
ete_tree.prune(keep_raw, preserve_branch_length=True)
print(f"  Pruned in {time.time()-t0:.1f}s")
print(f"  Leaves after pruning: {len(ete_tree.get_leaves())}")

PRUNED_TREE = f"{BASE}/data/bac120_r226_uniprot_pruned.tree"
ete_tree.write(format=1, outfile=PRUNED_TREE)
print(f"  Saved pruned tree: {PRUNED_TREE}")

tree = Phylo.read(PRUNED_TREE, "newick")


# -- 5. get leaf order from tree -------------------------------------------
print("\n=== 5. Extracting leaf order ===")
leaf_order_raw  = [cl.name for cl in tree.get_terminals()]
leaf_order      = [norm_id(n) for n in leaf_order_raw]
if len(set(leaf_order)) != len(leaf_order):
    seen = set()
    unique_raw = []
    unique_norm = []
    for raw, norm in zip(leaf_order_raw, leaf_order):
        if norm not in seen:
            seen.add(norm)
            unique_raw.append(raw)
            unique_norm.append(norm)
    leaf_order_raw = unique_raw
    leaf_order = unique_norm
print(f"  Leaf order length: {len(leaf_order)}")

leaf_set = set(leaf_order)
genome_list = [g for g in leaf_order
               if g in (set(hmm.index) & set(ipr.index)
                        & set(pub.index) & set(gtdb.index))]
print(f"  Final genome count: {len(genome_list)}")


# -- 6. architecture -> domain presence ------------------------------------
print("\n=== 6. Parsing domain presence ===")

def parse_arch(arch: str):
    """Return (ku, lig, pol, pe) booleans from Architecture string."""
    a = str(arch).lower()
    ku  = int("yes-nhej" in a or a == "ku-only")
    lig = int(a not in ("no-ligd-domains", "ku-only")
              and "no-ligd-domains" not in a)
    pol = int(("pol" in a and "no-pol" not in a) or "all" in a)
    pe  = int("all" in a)
    return ku, lig, pol, pe

tests = {
    "No-LigD-domains":          (0,0,0,0),
    "Ku-only":                   (1,0,0,0),
    "LigD-only-lig-no-pol":      (0,1,0,0),
    "LigD-only-lig-pol-diff-prot":(0,1,1,0),
    "LigD-only-all-diff-prot":   (0,1,1,1),
    "LigD-only-all-same-prot":   (0,1,1,1),
    "Yes-NHEJ-lig-no-pol":       (1,1,0,0),
    "Yes-NHEJ-lig-pol-diff-prot":(1,1,1,0),
    "Yes-NHEJ-all-diff-prot":    (1,1,1,1),
    "Yes-NHEJ-all-same-prot":    (1,1,1,1),
}
for arch, expected in tests.items():
    got = parse_arch(arch)
    status = "OK" if got == expected else f"FAIL (expected {expected})"
    print(f"  {arch}: {got} {status}")


# -- 7. build binary matrix ------------------------------------------------
print("\n=== 7. Building 12-column binary matrix ===")
cols = ["KU_PUB","KU_HMM","KU_IPR",
        "LIG_PUB","LIG_HMM","LIG_IPR",
        "POL_PUB","POL_HMM","POL_IPR",
        "PE_PUB","PE_HMM","PE_IPR"]
mat = pd.DataFrame(0, index=genome_list, columns=cols)

for gid in genome_list:
    ku, lig, pol, pe = parse_arch(hmm.loc[gid, "Architecture"])
    mat.loc[gid, "KU_HMM"]  = ku
    mat.loc[gid, "LIG_HMM"] = lig
    mat.loc[gid, "POL_HMM"] = pol
    mat.loc[gid, "PE_HMM"]  = pe

for gid in genome_list:
    ku, lig, pol, pe = parse_arch(ipr.loc[gid, "Architecture"])
    mat.loc[gid, "KU_IPR"]  = ku
    mat.loc[gid, "LIG_IPR"] = lig
    mat.loc[gid, "POL_IPR"] = pol
    mat.loc[gid, "PE_IPR"]  = pe

for gid in genome_list:
    row = pub.loc[gid]
    ku_id  = str(row["Ku_IDs"]).strip()
    lig_id = str(row["LigD_IDs"]).strip()
    ku  = int(ku_id.upper() not in ("NA-KU", "NA", "NAN", ""))
    lig = int(lig_id.upper() not in ("NA-LIGD", "NA", "NAN", ""))
    mat.loc[gid, "KU_PUB"]  = ku
    mat.loc[gid, "LIG_PUB"] = lig
    mat.loc[gid, "POL_PUB"] = lig
    mat.loc[gid, "PE_PUB"]  = lig

leaf_order_final = [g for g in leaf_order if g in mat.index]
mat = mat.loc[leaf_order_final]
leaf_order = leaf_order_final

print("  Column totals:")
for col in cols:
    print(f"    {col}: {mat[col].sum()}")

mat.to_csv(OUT_CSV)
print(f"  Saved: {OUT_CSV}")


# -- 8. parse GTDB taxonomy (Phylum + Class) -------------------------------
print("\n=== 8. Parsing GTDB taxonomy (Phylum/Class) ===")

def parse_taxon(tax_str: str):
    """Return (phylum, class) extracted from GTDB taxonomy string."""
    s = str(tax_str)
    phy = "Unknown"
    cls = "Unknown"
    for chunk in s.split(";"):
        chunk = chunk.strip()
        if chunk.startswith("p__"):
            v = chunk[3:].strip()
            phy = v if v else "Unknown"
        elif chunk.startswith("c__"):
            v = chunk[3:].strip()
            cls = v if v else "Unknown"
    return phy, cls

tax_rows = []
for gid in leaf_order:
    p, c = parse_taxon(gtdb.loc[gid, "gtdb_taxonomy"])
    tax_rows.append((gid, p, c))
tax_df = pd.DataFrame(tax_rows, columns=["gid", "Phylum", "Class"]).set_index("gid")
tax_df.to_csv(OUT_TAX)
print(f"  Unique Phyla : {tax_df['Phylum'].nunique()}")
print(f"  Unique Classes: {tax_df['Class'].nunique()}")
print("  Top 5 phyla:")
for p, n in tax_df["Phylum"].value_counts().head(5).items():
    print(f"    {p}: {n}")
print("  Top 5 classes:")
for c, n in tax_df["Class"].value_counts().head(5).items():
    print(f"    {c}: {n}")
print(f"  Saved taxonomy table: {OUT_TAX}")


def build_categorical_palette(categories, base_cmap_name="tab20",
                              extra_cmap_name="tab20b"):
    """Stable color mapping for a list of categories.

    Combines tab20 + tab20b + tab20c to provide up to 60 distinct
    colors.  Categories are sorted by frequency (most common first)
    so the visually dominant taxa get the most distinct colors.
    Always assigns a light grey to "Unknown".
    """
    pool = []
    for name in (base_cmap_name, "tab20b", "tab20c"):
        cmap = plt.get_cmap(name)
        for i in range(cmap.N):
            pool.append(to_hex(cmap(i)))
    palette = {}
    idx = 0
    for cat in categories:
        if cat == "Unknown":
            palette[cat] = "#D9D9D9"
            continue
        palette[cat] = pool[idx % len(pool)]
        idx += 1
    return palette

# order categories by frequency so the largest groups get visually
# distinct colors first
phyla_ordered   = list(tax_df["Phylum"].value_counts().index)
classes_ordered = list(tax_df["Class"].value_counts().index)
phylum_palette  = build_categorical_palette(phyla_ordered)
class_palette   = build_categorical_palette(classes_ordered)


# -- 9. build figure -------------------------------------------------------
print("\n=== 9. Building figure ===")

N = len(leaf_order)
FIG_W   = 16
FIG_H   = max(10, N * 0.06)

# horizontal layout fractions
TREE_FRAC      = 0.30
LABEL_FRAC     = 0.09
HEAT_FRAC      = 0.45      # 12 NHEJ tracks
GAP_TAX        = 0.012     # gap between NHEJ and taxonomy block
TAX_FRAC       = 0.05      # 2 taxonomy tracks
TRACK_W        = HEAT_FRAC / 12
TAX_TRACK_W    = TAX_FRAC / 2

# colours per domain (3 shades: dark=PUB, mid=HMM, light=IPR)
DOMAIN_COLORS = {
    "KU":  ("#1565C0", "#1E88E5", "#90CAF9"),
    "LIG": ("#2E7D32", "#43A047", "#A5D6A7"),
    "POL": ("#B71C1C", "#E53935", "#EF9A9A"),
    "PE":  ("#E65100", "#FB8C00", "#FFCC80"),
}

fig = plt.figure(figsize=(FIG_W, FIG_H))
fig.patch.set_facecolor("white")

left_tree   = 0.02
left_label  = left_tree + TREE_FRAC + 0.01
left_tracks = left_label + LABEL_FRAC + 0.005
bottom      = 0.04
height_ax   = 0.92

ax_tree  = fig.add_axes([left_tree, bottom, TREE_FRAC, height_ax])
ax_label = fig.add_axes([left_label, bottom, LABEL_FRAC, height_ax])

track_axes = {}
x_cursor = left_tracks
domain_order = ["KU","LIG","POL","PE"]
source_order = ["PUB","HMM","IPR"]
col_order = [f"{d}_{s}" for d in domain_order for s in source_order]

for col in col_order:
    ax = fig.add_axes([x_cursor, bottom, TRACK_W * 0.90, height_ax])
    track_axes[col] = ax
    x_cursor += TRACK_W

# add taxonomy tracks (Phylum, Class) after a gap
x_cursor += GAP_TAX
tax_axes = {}
for tname in ("Phylum", "Class"):
    ax = fig.add_axes([x_cursor, bottom, TAX_TRACK_W * 0.90, height_ax])
    tax_axes[tname] = ax
    x_cursor += TAX_TRACK_W


# -- draw tree ------------------------------------------------------------
print(f"  Drawing tree ({N} leaves)...")
ax_tree.set_xlim(-0.02, 1.05)
ax_tree.set_ylim(0, N)
ax_tree.axis("off")

y_pos = {leaf_order[i]: i + 0.5 for i in range(N)}

def get_max_depth(clade):
    if clade.is_terminal():
        return 0.0
    return max(
        (c.branch_length or 0.0) + get_max_depth(c)
        for c in clade.clades
    )

def draw_subtree(clade, x, ax, x_scale, y_pos, color="#333333"):
    bl = clade.branch_length if clade.branch_length else 0.0
    x_new = x + bl * x_scale
    if clade.is_terminal():
        return x_new, y_pos[norm_id(clade.name)]
    child_x_ends, child_y_ends = [], []
    for child in clade.clades:
        cx, cy = draw_subtree(child, x_new, ax, x_scale, y_pos, color)
        child_x_ends.append(cx)
        child_y_ends.append(cy)
    y_lo = min(child_y_ends)
    y_hi = max(child_y_ends)
    ax.plot([x_new, x_new], [y_lo, y_hi],
            color=color, lw=0.4, solid_capstyle="round")
    y_mid = (y_lo + y_hi) / 2.0
    for cx, cy in zip(child_x_ends, child_y_ends):
        ax.plot([x_new, cx], [cy, cy],
                color=color, lw=0.4, solid_capstyle="round")
    return x_new, y_mid

root = tree.root
max_depth = get_max_depth(root) or 1.0
x_scale = 0.95 / max_depth

t_draw = time.time()
draw_subtree(root, 0, ax_tree, x_scale, y_pos)
print(f"  Tree drawn in {time.time()-t_draw:.1f}s")
ax_tree.set_yticks([])
ax_tree.set_xticks([])


# -- draw genome labels ---------------------------------------------------
ax_label.set_xlim(0, 1)
ax_label.set_ylim(0, N)
ax_label.axis("off")
if N <= 1500:
    fs = max(2.0, min(5.0, 300.0 / N))
    for i, gid in enumerate(leaf_order):
        ax_label.text(0.0, i + 0.5, gid, va="center", ha="left",
                      fontsize=fs, fontfamily="monospace",
                      color="#222222")
else:
    ax_label.text(0.5, N/2, f"({N} genomes)", va="center", ha="center",
                  fontsize=7, color="#666666")


# -- draw 12 binary tracks ------------------------------------------------
print("  Drawing 12 binary tracks...")
for col in col_order:
    domain, source = col.split("_")
    src_idx = {"PUB": 0, "HMM": 1, "IPR": 2}[source]
    color = DOMAIN_COLORS[domain][src_idx]

    ax = track_axes[col]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, N)
    ax.axis("off")
    vals = mat[col].values
    for i, v in enumerate(vals):
        if v == 1:
            ax.add_patch(mpatches.Rectangle(
                (0.05, i + 0.1), 0.90, 0.80,
                color=color, lw=0, zorder=2
            ))
    ax.text(0.5, N + N*0.01, col.replace("_","\n"),
            ha="center", va="bottom", fontsize=4.5,
            fontweight="bold", color=color,
            transform=ax.transData)

# domain-group separator lines (between KU|LIG, LIG|POL, POL|PE)
for di in range(1, 4):
    col_at_boundary = col_order[di * 3]
    ax_at = track_axes[col_at_boundary]
    ax_at.axvline(0, color="black", lw=0.8, alpha=0.4, zorder=5)


# -- draw 2 taxonomy tracks (Phylum, Class) -------------------------------
print("  Drawing 2 taxonomy tracks (Phylum, Class)...")
phy_vals = tax_df.loc[leaf_order, "Phylum"].tolist()
cls_vals = tax_df.loc[leaf_order, "Class"].tolist()
tax_data = {"Phylum": (phy_vals, phylum_palette),
            "Class":  (cls_vals, class_palette)}

for tname, (vals, palette) in tax_data.items():
    ax = tax_axes[tname]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, N)
    ax.axis("off")
    for i, v in enumerate(vals):
        c = palette.get(v, "#D9D9D9")
        ax.add_patch(mpatches.Rectangle(
            (0.05, i), 0.90, 1.0,
            color=c, lw=0, zorder=2
        ))
    ax.text(0.5, N + N*0.01, tname,
            ha="center", va="bottom", fontsize=5.5,
            fontweight="bold", color="#222222",
            transform=ax.transData)
    ax.axvline(0, color="black", lw=0.8, alpha=0.4, zorder=5)


# -- legends --------------------------------------------------------------
# Domain x source legend
domain_handles = []
for domain, (c0, c1, c2) in DOMAIN_COLORS.items():
    domain_handles += [
        mpatches.Patch(facecolor=c0, label=f"{domain} PUB"),
        mpatches.Patch(facecolor=c1, label=f"{domain} HMM"),
        mpatches.Patch(facecolor=c2, label=f"{domain} IPR"),
    ]
leg1 = fig.legend(handles=domain_handles,
                  loc="lower center",
                  ncol=4, fontsize=6, framealpha=0.9,
                  bbox_to_anchor=(0.30, 0.005),
                  title="NHEJ domain x source", title_fontsize=7)

# Phylum legend (top phyla; group rare ones into "Other")
def make_taxon_handles(values, palette, max_show=18, label="Phylum"):
    counts = pd.Series(values).value_counts()
    handles = []
    shown = 0
    other_n = 0
    for cat, n in counts.items():
        if shown < max_show:
            handles.append(mpatches.Patch(
                facecolor=palette.get(cat, "#D9D9D9"),
                label=f"{cat} (n={n})"))
            shown += 1
        else:
            other_n += n
    if other_n > 0:
        handles.append(mpatches.Patch(
            facecolor="#999999",
            label=f"Other {label.lower()} (n={other_n})"))
    return handles

phy_handles = make_taxon_handles(phy_vals, phylum_palette,
                                 max_show=18, label="Phylum")
cls_handles = make_taxon_handles(cls_vals, class_palette,
                                 max_show=18, label="Class")

leg2 = fig.legend(handles=phy_handles,
                  loc="lower center", ncol=3, fontsize=5,
                  framealpha=0.9,
                  bbox_to_anchor=(0.66, 0.005),
                  title="Phylum (top groups)", title_fontsize=6)

leg3 = fig.legend(handles=cls_handles,
                  loc="lower center", ncol=3, fontsize=5,
                  framealpha=0.9,
                  bbox_to_anchor=(0.88, 0.005),
                  title="Class (top groups)", title_fontsize=6)

fig.add_artist(leg1)
fig.add_artist(leg2)

# title
fig.suptitle(
    f"Phylogenetic distribution of NHEJ domains across {N} bacterial genomes\n"
    f"(GTDB bac120 R226 tree pruned to intersection of HMM, InterPro, "
    f"published data, and GTDB taxonomy)",
    fontsize=9, y=0.995, ha="center", va="top"
)

print(f"\n=== 10. Saving figure ({FIG_W}x{FIG_H:.1f} in) ===")
plt.savefig(OUT_PDF, bbox_inches="tight", dpi=300)
plt.savefig(OUT_PNG, bbox_inches="tight", dpi=150)
plt.close()
print(f"  Saved PDF: {OUT_PDF}")
print(f"  Saved PNG: {OUT_PNG}")

print("\n=== DONE ===")
