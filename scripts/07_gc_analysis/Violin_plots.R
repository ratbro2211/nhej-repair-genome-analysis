setwd("/Users/surajrathi/Downloads/Repair_enzyme/pfam/ku")
# -------------------------------
# 1. LOAD LIBRARIES
# -------------------------------
library(ggplot2)
library(dplyr)

# -------------------------------
# 2. LOAD DATA
# -------------------------------
df <- read.delim(
  "tree_leaf_PF02735_presence_with_GC_FINAL.tsv",
  stringsAsFactors = FALSE
)

# -------------------------------
# 3. NORMALIZE COLUMNS
# -------------------------------

# ku_present → character
df$ku_present <- as.character(df$ku_present)

# Convert to readable labels
df$Ku_status <- ifelse(
  tolower(df$ku_present) == "true",
  "Ku",
  "No Ku"
)

# GC to numeric (not found → NA automatically)
df$gc_percentage <- suppressWarnings(as.numeric(df$gc_percentage))

# -------------------------------
# 4. SANITY CHECK (IMPORTANT)
# -------------------------------
cat("\nCounts (ALL rows):\n")

print(table(df$Ku_status))

cat("\nCounts with usable GC:\n")
print(table(df$Ku_status, !is.na(df$gc_percentage)))

# -------------------------------
# 5. VIOLIN PLOT
# -------------------------------
p <- ggplot(
  df,
  aes(x = Ku_status, y = gc_percentage)
) +
  geom_violin(
    fill = "grey80",
    color = "black",
    na.rm = TRUE,
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.08,
    fill = "black",
    outlier.shape = NA,
    na.rm = TRUE
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3,
    fill = "white",
    na.rm = TRUE
  ) +
  labs(
    x = "",
    y = "GC content (%)",
    title = "GC Content Distribution by Ku Presence"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# -------------------------------
# 6. SAVE
# -------------------------------
pdf("Results/GC_content_Ku_violin.pdf", width = 5, height = 6)
print(p)
dev.off()

cat("\n✅ Violin plot saved: Results/GC_content_Ku_violin.pdf\n")



library(dplyr)

# -------------------------------
# 1. LOAD DATA
# -------------------------------
df <- read.delim(
  "tree_leaf_PF02735_presence_with_GC_FINAL.tsv",
  stringsAsFactors = FALSE
)

# Normalize Ku status
df$Ku_status <- ifelse(
  tolower(df$ku_present) == "true",
  "Ku",
  "No Ku"
)

# Convert GC safely
df$gc_percentage_num <- suppressWarnings(as.numeric(df$gc_percentage))


# -------------------------------
# 2. OVERALL DATASET SUMMARY
# -------------------------------
overall_summary <- data.frame(
  Total_genomes = nrow(df),
  Ku_positive = sum(df$Ku_status == "Ku"),
  Ku_negative = sum(df$Ku_status == "No Ku"),
  stringsAsFactors = FALSE
)

print(overall_summary)


# -------------------------------
# 3. GC METADATA AVAILABILITY
# -------------------------------
gc_summary <- data.frame(
  Total_genomes = nrow(df),
  GC_available = sum(!is.na(df$gc_percentage_num)),
  GC_missing = sum(is.na(df$gc_percentage_num)),
  stringsAsFactors = FALSE
)

print(gc_summary)


# -------------------------------
# 4. GC AVAILABILITY BY Ku STATUS
# -------------------------------
gc_by_ku <- df %>%
  mutate(GC_status = ifelse(is.na(gc_percentage_num), "Missing", "Available")) %>%
  count(Ku_status, GC_status) %>%
  tidyr::pivot_wider(
    names_from = GC_status,
    values_from = n,
    values_fill = 0
  )

print(gc_by_ku)


# -------------------------------
# 5. GC STATISTICS (AVAILABLE ONLY)
# -------------------------------
gc_stats <- df %>%
  filter(!is.na(gc_percentage_num)) %>%
  group_by(Ku_status) %>%
  summarise(
    N = n(),
    Mean_GC = mean(gc_percentage_num),
    Median_GC = median(gc_percentage_num),
    SD_GC = sd(gc_percentage_num),
    Min_GC = min(gc_percentage_num),
    Max_GC = max(gc_percentage_num)
  )

print(gc_stats)


# -------------------------------
# 6. SAVE ALL SUMMARIES
# -------------------------------
write.csv(
  overall_summary,
  "Results/summary_overall_dataset.csv",
  row.names = FALSE
)

write.csv(
  gc_summary,
  "Results/summary_gc_availability.csv",
  row.names = FALSE
)

write.csv(
  gc_by_ku,
  "Results/summary_gc_by_Ku.csv",
  row.names = FALSE
)

write.csv(
  gc_stats,
  "Results/summary_gc_statistics_by_Ku.csv",
  row.names = FALSE
)

cat("\n✅ All summary tables saved in Results/\n")

