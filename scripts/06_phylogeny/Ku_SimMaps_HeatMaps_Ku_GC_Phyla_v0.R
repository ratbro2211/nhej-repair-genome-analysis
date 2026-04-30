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

library(plotrix)

# ---- 1. Posterior probabilities ----
pp <- summary(simmap)$ace[1:tree$Nnode, ]

# ---- 2. Density map ----
dmap <- densityMap(simmap, plot=FALSE)

# ---- 3. Colors ----
cols <- setNames(c("blue","red"), c("0","1"))

# ---- 4. PDF ----
pdf("ku_gc_phylogeny_heatmap.pdf", width=16, height=48)
par(mar=c(5,1,1,30))

# ---- 5. Plot tree ----
plot(dmap,
     type="phylogram",
     ftype="off",
     lwd=0.6,
     lend=1,
     legend=FALSE,
     xlim=c(0, max(nodeHeights(tree))*2),
     ylim=c(0, Ntip(tree)*2))

# ---- 6. Coordinates ----
obj <- get("last_plot.phylo", envir=.PlotPhyloEnv)
n <- Ntip(tree)
h <- max(nodeHeights(tree))

# ---- 7. Filter pies (correct logic) ----
keep <- which(apply(pp, 1, function(x){
  m <- max(x)
  (m > 0.8) | (m < 0.1)
}))

ii <- keep[order(rowSums(pp[keep,]^2), decreasing=TRUE)]

# ---- 8. Draw pies ----
par(xpd=TRUE)

for(i in ii){
  plotrix::floating.pie(
    xpos = obj$xx[i+n],
    ypos = obj$yy[i+n],
    x = pp[i,],
    radius = 0.0018 * h,
    col = cols,
    border = NA
  )
}

# =========================
# DEFINE POSITIONS
# =========================

x_tree_end <- max(obj$xx)

# ---- LABEL POSITIONS (BETWEEN TREE & HEATMAPS) ----
x_label1 <- x_tree_end + 0.03   # GCF
x_label2 <- x_tree_end + 0.18   # taxonomy

# ---- HEATMAP START ----
x_offset <- x_tree_end + 0.8

# =========================
# LABELS FIRST (CRITICAL)
# =========================

show_idx <- seq(1, length(obj$yy), by=20)

text(x_label1, obj$yy[show_idx],
     df$tree_leaf[show_idx],
     cex=0.25, pos=4)

text(x_label2, obj$yy[show_idx],
     df$genus_species[show_idx],
     cex=0.25, pos=4)

# =========================
# HEATMAPS
# =========================

# Ku heatmap
for(i in 1:length(obj$yy)){
  rect(x_offset,
       obj$yy[i]-0.5,
       x_offset+0.04,
       obj$yy[i]+0.5,
       col=ifelse(ku[i]==1,"red","white"),
       border=NA)
}

# GC heatmap
gc_vec <- as.numeric(as.character(df$gc_percentage))
gc_vec[is.na(gc_vec)] <- median(gc_vec, na.rm=TRUE)

gc_cols <- colorRampPalette(c("blue","yellow","red"))(100)
gc_scaled <- as.numeric(cut(gc_vec, breaks=100))

for(i in 1:length(obj$yy)){
  rect(x_offset+0.05,
       obj$yy[i]-0.5,
       x_offset+0.09,
       obj$yy[i]+0.5,
       col=gc_cols[gc_scaled[i]],
       border=NA)
}

# Phylum heatmap
df$phylum <- sub(".*;p__", "p__", df$gtdb_taxonomy)
df$phylum <- sub(";.*", "", df$phylum)

phy_levels <- unique(df$phylum)
phy_cols <- setNames(rainbow(length(phy_levels)), phy_levels)

for(i in 1:length(obj$yy)){
  rect(x_offset+0.10,
       obj$yy[i]-0.5,
       x_offset+0.14,
       obj$yy[i]+0.5,
       col=phy_cols[df$phylum[i]],
       border=NA)
}

# =========================
# COLUMN LABELS (ABOVE HEATMAPS)
# =========================

y_top <- max(obj$yy)
label_y <- y_top + 1

text(x_offset+0.02, label_y, "Ku",    cex=0.8, font=2)
text(x_offset+0.07, label_y, "GC%",   cex=0.8, font=2)
text(x_offset+0.12, label_y, "Phyla", cex=0.8, font=2)

# =========================
# LEGENDS (RIGHT SIDE ONLY)
# =========================

x_leg <- x_offset + 0.25

# GC legend
gc_leg_cols <- colorRampPalette(c("blue","yellow","red"))(5)

legend(x = x_leg,
       y = y_top * 0.95,
       legend = round(seq(min(gc_vec), max(gc_vec), length.out=5),1),
       fill = gc_leg_cols,
       title = "GC %",
       border = NA,
       bty = "n",
       cex = 0.7)

# Ku legend
legend(x = x_leg,
       y = y_top * 0.75,
       legend = c("Ku absent", "Ku present"),
       fill = c("blue","red"),
       border = NA,
       bty = "n",
       cex = 0.7)

# Gain/loss text
text(x = x_leg,
     y = y_top * 0.60,
     labels = "Blue = absence\nRed = presence",
     cex = 0.6,
     adj = 0)

# Phylum legend (top 10)
top_phy <- names(sort(table(df$phylum), decreasing=TRUE))[1:10]

legend(x = x_leg,
       y = min(obj$yy) + (y_top - min(obj$yy))*0.15,
       legend = top_phy,
       fill = phy_cols[top_phy],
       border = NA,
       bty = "n",
       cex = 0.55,
       title = "Top Phyla")

dev.off()

