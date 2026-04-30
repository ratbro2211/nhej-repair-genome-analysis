# ================================
# SIMMAP (EVENT DENSITY) + HEATMAPS + CUSTOM ROW LABELS (PDF)
# ================================

library(phytools)
library(ape)
library(data.table)

# ---- 0. File paths ----
BASE = "/Users/surajrathi/Downloads/Repair_enzyme/pfam/ku"
TREE_FILE <- file.path(BASE, "bac120_r226_uniprot_pruned.tree")
GC_FILE   <- file.path(BASE, "tree_leaf_PF02735_presence_with_GC_FINAL.tsv")
GTDB_FILE <- file.path(BASE, "uniprot_proteomes_with_GTDB_R226.tsv")

# ---- 1. Load data ----
tree <- read.tree(TREE_FILE)

gc <- read.table(GC_FILE, sep="\t", header=TRUE, stringsAsFactors=FALSE)
gtdb <- fread(GTDB_FILE, sep="\t", header=TRUE, data.table=FALSE)

# ---- 2. Clean column names ----
colnames(gc)   <- tolower(trimws(colnames(gc)))
colnames(gtdb) <- tolower(trimws(colnames(gtdb)))

# ---- 3. Normalize IDs ----
normalize_id <- function(x){
  x <- trimws(x)
  x <- sub("^RS_|^GB_", "", x)
  x <- sub("^GCF_|^GCA_", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  return(x)
}

gc$tree_leaf   <- normalize_id(gc$tree_leaf)
gtdb$accession <- normalize_id(gtdb$accession)
tree$tip.label <- normalize_id(tree$tip.label)

# ---- 4. Remove blanks ----
gc   <- gc[!is.na(gc$tree_leaf)   & gc$tree_leaf   != "", ]
gtdb <- gtdb[!is.na(gtdb$accession) & gtdb$accession != "", ]

# ---- 5. Merge ----
df <- merge(gc, gtdb, by.x="tree_leaf", by.y="accession")

# ---- 6. Prepare Ku trait ----
ku <- setNames(df$ku_present, df$tree_leaf)
ku <- ifelse(ku=="True", 1, 0)

# ---- 7. Match with tree ----
common <- intersect(tree$tip.label, names(ku))
tree <- drop.tip(tree, setdiff(tree$tip.label, common))
ku   <- ku[common]

valid <- !is.na(ku)
tree <- drop.tip(tree, setdiff(tree$tip.label, names(ku)[valid]))
ku   <- ku[valid]

ku <- factor(ku, levels=c(0,1))

# ---- 8. Reorder df ----
df <- df[df$tree_leaf %in% tree$tip.label, ]
df <- df[match(tree$tip.label, df$tree_leaf), ]

# ---- 9. Extract genus_species ----
df$genus_species <- sub(".*;g__", "", df$gtdb_taxonomy)
df$genus_species <- sub(";s__", "_", df$genus_species)
df$genus_species <- sub(";.*", "", df$genus_species)

# ---- 10. SIMMAP ----
Q <- matrix(c(-1,1,1,-1),2,2)
rownames(Q) <- colnames(Q) <- c("0","1")

set.seed(1)
simmap <- make.simmap(tree, ku, model=Q, nsim=100)

# ================================
# CORRECT phytools-style: densityMap + pies + heatmaps
# ================================

library(phytools)
library(ape)
library(data.table)
library(plotrix)

# ---- SIMMAP already computed ----
# simmap

# ---- 1. Posterior probs ----
pp <- summary(simmap)$ace[1:tree$Nnode, ]

# ---- 2. Density map (for coordinates ONLY) ----
dmap <- densityMap(simmap, plot=FALSE)

# ---- 3. Colors ----
# same gradient for branches AND pies
cols <- setNames(c("blue","red"), c("0","1"))

# ---- 4. PDF ----
pdf("ku_phytools_style.pdf", width=16, height=48)
par(mar=c(5,1,1,25))

# ---- 5. Plot density map (RECTANGULAR, no labels) ----
plot(dmap,
     type="phylogram",   # rectangular (NOT fan)
     ftype="off",
     lwd=1,
     legend=FALSE,
     xlim=c(0, max(nodeHeights(tree))*2))

# ---- 6. Get coordinates ----
obj <- get("last_plot.phylo", envir=.PlotPhyloEnv)
n <- Ntip(tree)
h <- max(nodeHeights(tree))

# ---- 7. Draw node pies (same logic as example) ----
ii <- order(rowSums(pp^2), decreasing=TRUE)

par(xpd=TRUE)

for(i in ii){
  plotrix::floating.pie(
    xpos = obj$xx[i+n],
    ypos = obj$yy[i+n],
    x = pp[i,],
    radius = 0.003 * h,   # << smaller pies
    col = cols,           # << SAME colors as branches
    border = NA
  )
}

# ---- 8. Heatmaps ----
x_offset <- max(obj$xx) * 1.1

# Ku
for(i in 1:length(obj$yy)){
  rect(x_offset,
       obj$yy[i]-0.15,
       x_offset+0.04,
       obj$yy[i]+0.15,
       col=ifelse(ku[i]==1,"red","white"),
       border=NA)
}

# GC
gc_vec <- as.numeric(as.character(df$gc_percentage))
gc_vec[is.na(gc_vec)] <- median(gc_vec, na.rm=TRUE)

gc_cols <- colorRampPalette(c("blue","yellow","red"))(100)
gc_scaled <- as.numeric(cut(gc_vec, breaks=100))

for(i in 1:length(obj$yy)){
  rect(x_offset+0.05,
       obj$yy[i]-0.15,
       x_offset+0.09,
       obj$yy[i]+0.15,
       col=gc_cols[gc_scaled[i]],
       border=NA)
}

# Phylum
df$phylum <- sub(".*;p__", "p__", df$gtdb_taxonomy)
df$phylum <- sub(";.*", "", df$phylum)

phy_levels <- unique(df$phylum)
phy_cols <- setNames(rainbow(length(phy_levels)), phy_levels)

for(i in 1:length(obj$yy)){
  rect(x_offset+0.10,
       obj$yy[i]-0.15,
       x_offset+0.14,
       obj$yy[i]+0.15,
       col=phy_cols[df$phylum[i]],
       border=NA)
}

# ---- 9. Labels ----
df$genus_species <- sub(".*;g__", "", df$gtdb_taxonomy)
df$genus_species <- sub(";s__", "_", df$genus_species)
df$genus_species <- sub(";.*", "", df$genus_species)

text(x_offset+0.16, obj$yy, labels=df$tree_leaf, cex=0.25, pos=4)
text(x_offset+0.35, obj$yy, labels=df$genus_species, cex=0.25, pos=4)

# ---- 10. Legend ----
legend("topright",
       legend=c("Ku absent","Ku present"),
       fill=cols,
       border=NA,
       bty="n")

dev.off()