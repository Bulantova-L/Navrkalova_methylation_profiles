library(biomaRt)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

setwd("z:/share/710000-CEITEC/713000-cmm/713016-bioit/base/workspace/lucka/Navrkalova_cfDNA/methylation_profiles")

read_cns_file <- function(file_path) {
  sample_id <- sub("_cns_mapped\\.bed$", "", basename(file_path))
  sample_id <- sub("_methylation$", "", sample_id)

  dt <- fread(file_path, header = FALSE)

  if (ncol(dt) != 6) {
    stop("Unexpected number of columns in: ", file_path, " (got ", ncol(dt), ")")
  }

  out <- data.table(
    sample = sample_id,
    chr = dt[[1]],
    pos = dt[[2]],
    end = dt[[3]],
    gene = dt[[4]],
    methylation = as.numeric(dt[[5]]),
    coverage = as.numeric(dt[[6]])
  )

  out[, probe_id := paste0(chr, ":", pos)]

  out
}

files <- list.files(
  "results",
  pattern = "_cns_mapped\\.bed$",
  recursive = TRUE,
  full.names = TRUE
)

cns_long <- rbindlist(lapply(files, read_cns_file), fill = TRUE)

gene_summary <- cns_long[
  !is.na(methylation),
  .(
    mean_methylation = mean(methylation, na.rm = TRUE),
    median_methylation = median(methylation, na.rm = TRUE),
    n_cpg = .N
  ),
  by = .(sample, gene)
]

gene_summary[order(sample, gene)]

gene_wide <- dcast(
  gene_summary,
  sample ~ gene,
  value.var = "mean_methylation"
)

# This is the main table for:heatmaps, PCA, clustering, risk scoring

# 4) Compute a simple CNS methylation score

cns_genes <- c(
  "CXCR4", "CXCL12", "S1PR1", "S1PR2",
  "MMP9", "CLDN5", "TJP1",
  "CD274", "HLA-A", "B2M",
  "MYD88", "CD79B", "PIM1", "PRDM1"
)

score_table <- gene_summary[gene %in% cns_genes]

# sample-level score: mean across all target genes
sample_score <- score_table[
  ,
  .(cns_mean_methylation = mean(mean_methylation, na.rm = TRUE)),
  by = sample
][order(-cns_mean_methylation)]

score_split <- gene_summary[gene %in% cns_genes][
  ,
  .(
    homing = mean(mean_methylation[gene %in% c("CXCR4", "CXCL12", "S1PR1", "S1PR2")], na.rm = TRUE),
    bbb = mean(mean_methylation[gene %in% c("MMP9", "CLDN5", "TJP1")], na.rm = TRUE),
    immune_escape = mean(mean_methylation[gene %in% c("CD274", "HLA-A", "B2M")], na.rm = TRUE),
    oncogenic = mean(mean_methylation[gene %in% c("MYD88", "CD79B", "PIM1", "PRDM1")], na.rm = TRUE)
  ),
  by = sample
]
# CNS homing, blood-brain barrier, immune escape, and oncogenic methylation scores

score_split

# 5) Plot the results
# Heatmap

library(ggplot2)

heatmap_plot <- ggplot(score_table, aes(x = gene, y = sample, fill = mean_methylation)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL, fill = "Meth", title = "CNS gene methylation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sample_score_plot <- ggplot(sample_score, aes(x = reorder(sample, cns_mean_methylation), y = cns_mean_methylation)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Mean CNS methylation", title = "Sample-level CNS methylation score") +
  theme_minimal()

ggsave("analyses/03_cns_genes/plots/gene_methylation_heatmap.png", plot = heatmap_plot, width = 10, height = 6, dpi = 300)
ggsave("analyses/03_cns_genes/plots/gene_methylation_sample_score.png", plot = sample_score_plot, width = 8, height = 5, dpi = 300)


