## Analysis of methylated CpG calls from modkit pileups:
##  1. bar chart of methylated CpG counts per sample
##  2. Venn/Euler diagram of shared methylated CpG positions
##     (cfDNA_VN = healthy, the rest = patient samples)
##  3. PCA on genome-wide 100kb methylation windows
##     - "absolute": window mean of per-CpG %-modified (mean_pct, as before)
##     - "relative": methylated CpGs / covered CpGs per window, which
##       normalizes out cross-sample differences in overall coverage
##
## Expects the outputs of 01_count_methylated_cpgs.sh in `analyses/01_cpg_counts_windows/output/`:
##   output/methylated_cpg_counts.tsv
##   output/{sample}_meth_cpgs.bed   (chr, start, end)
##   output/{sample}_windows.bed     (chr, start, end, mean_pct, n_cpgs)

setwd("/storage/713000-cmm/713016-bioit/base/workspace/lucka/Navrkalova_cfDNA/methylation_profiles")

library(data.table)
library(ggplot2)
library(polyclip)
library(eulerr)
library(polylabelr)

# ---- config -----------------------------------------------------------------
derived_dir <- "analyses/01_cpg_counts_windows/output"
plots_dir   <- "analyses/01_cpg_counts_windows/plots"
samples     <- c("cfDNA_BAR", "cfDNA_LOS", "cfDNA_MUD", "cfDNA_Pilot1", "cfDNA_VN")
healthy_sample <- "cfDNA_VN"
min_cpgs_per_window <- 5   # drop windows with too few covered CpGs to trust the mean
window_size <- 100000      # must match WINDOW_SIZE in 01b_bedtools_windows.sh
# -------------------------------------------------------------------------------

dir.create(plots_dir, showWarnings = FALSE)

group_of <- function(s) ifelse(s == healthy_sample, "Healthy", "Patient")

## 1. Methylated CpG counts per sample --------------------------------------
counts <- fread(file.path(derived_dir, "methylated_cpg_counts.tsv"))
counts[, group := group_of(sample)]
print(counts)

p_counts <- ggplot(counts, aes(x = reorder(sample, -n_methylated_cpgs), y = n_methylated_cpgs, fill = group)) +
  geom_col() +
  labs(x = NULL, y = "# methylated CpGs", title = "Methylated CpGs per sample") +
  theme_minimal()
ggsave(file.path(plots_dir, "methylated_cpg_counts.png"), p_counts, width = 6, height = 4)

## 2. Venn / Euler diagram of shared methylated CpG positions ----------------
chrom_levels <- paste0("chr", c(1:22, "X", "Y", "M"))

read_meth_keys <- function(s) {
  f <- file.path(derived_dir, paste0(s, "_meth_cpgs.bed"))
  d <- fread(f, select = 1:2, col.names = c("chr", "start"))
  chr_code <- match(d$chr, chrom_levels)
  pos_id <- unique(as.numeric(chr_code) * 1e10 + d$start)
  data.table(pos_id = pos_id, sample = s, present = TRUE)
}

dt_long <- rbindlist(lapply(samples, read_meth_keys))
wide <- dcast(dt_long, pos_id ~ sample, value.var = "present", fill = FALSE)
for (s in samples) wide[[s]] <- as.logical(wide[[s]])

combo_counts <- wide[, .N, by = samples]
setorder(combo_counts, -N)
fwrite(combo_counts, file.path(derived_dir, "venn_combo_counts.tsv"), sep = "\t")
print(combo_counts)

euler_input <- setNames(numeric(0), character(0))
for (i in seq_len(nrow(combo_counts))) {
  present <- samples[unlist(combo_counts[i, ..samples])]
  if (length(present) == 0) next
  euler_input[paste(present, collapse = "&")] <- combo_counts$N[i]
}

if (requireNamespace("eulerr", quietly = TRUE)) {
  library(eulerr)
  fit <- euler(euler_input)
  png(file.path(plots_dir, "venn_methylated_cpgs.png"), width = 1400, height = 1400, res = 200)
  print(plot(fit, quantities = TRUE, main = "Shared methylated CpGs"))
  dev.off()
} else {
  message("Package 'eulerr' not installed (install.packages(\"eulerr\")) - ",
          "falling back to an UpSet-style bar plot of intersection sizes.")
  combo_df <- combo_counts
  combo_df$combo <- apply(combo_df[, ..samples], 1, function(r) paste(samples[r], collapse = "&"))
  combo_df <- combo_df[combo != ""]
  p_upset <- ggplot(combo_df, aes(x = reorder(combo, -N), y = N)) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = "# shared methylated CpGs",
         title = "Methylated CpG overlap across samples (UpSet-style)") +
    theme_minimal()
  ggsave(file.path(plots_dir, "venn_methylated_cpgs_upset.png"), p_upset, width = 7, height = 5)
}

## 3. PCA on genome-wide 100kb methylation windows ---------------------------

# Per-window count of methylated CpGs (absolute, from {sample}_meth_cpgs.bed),
# binned into the same fixed windows as {sample}_windows.bed. Needed to turn
# the absolute methylated-CpG count into a coverage-normalized fraction.
read_meth_window_counts <- function(s) {
  f <- file.path(derived_dir, paste0(s, "_meth_cpgs.bed"))
  d <- fread(f, select = 1:2, col.names = c("chr", "start"))
  d[, win_start := (start %/% window_size) * window_size]
  d[, win_id := paste(chr, win_start, sep = "_")]
  d[, .(n_methylated = .N), by = win_id]
}

read_windows <- function(s) {
  f <- file.path(derived_dir, paste0(s, "_windows.bed"))
  # bedtools map fills windows with no overlapping CpGs as "."
  d <- fread(f, col.names = c("chr", "start", "end", "mean_pct", "n_cpgs"), na.strings = ".")
  d[, win_id := paste(chr, start, sep = "_")]
  d <- d[!is.na(n_cpgs) & n_cpgs >= min_cpgs_per_window, .(win_id, mean_pct, n_cpgs)]

  meth <- read_meth_window_counts(s)
  d <- merge(d, meth, by = "win_id", all.x = TRUE)
  d[is.na(n_methylated), n_methylated := 0]
  # relative methylation: fraction of THIS sample's covered CpGs in the window
  # that were called methylated - normalizes away sample-to-sample coverage
  # differences that inflate absolute methylated-CpG counts in deeper samples
  d[, relative_pct := 100 * n_methylated / n_cpgs]
  d[, .(win_id, mean_pct, relative_pct)]
}

per_sample <- lapply(samples, function(s) read_windows(s))
names(per_sample) <- samples

build_matrix <- function(value_col) {
  cols <- lapply(samples, function(s) {
    d <- per_sample[[s]][, c("win_id", value_col), with = FALSE]
    setnames(d, value_col, s)
    d
  })
  merged <- Reduce(function(x, y) merge(x, y, by = "win_id", all = FALSE), cols)
  mat <- as.matrix(merged[, ..samples])
  rownames(mat) <- merged$win_id
  mat[complete.cases(mat), , drop = FALSE]
}

run_pca_plot <- function(mat, title, out_file) {
  cat(sprintf("%s: %d windows x %d samples (windows with >=%d CpGs covered in every sample)\n",
              title, nrow(mat), ncol(mat), min_cpgs_per_window))

  pca <- prcomp(t(mat), scale. = TRUE)
  var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  pca_df <- data.frame(
    sample = colnames(mat),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group = group_of(colnames(mat))
  )

  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
    geom_point(size = 4) +
    geom_text(aes(label = sample), vjust = -1, show.legend = FALSE, size = 3) +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    theme_minimal()
  ggsave(file.path(plots_dir, out_file), p_pca, width = 6, height = 5)
  invisible(pca)
}

# absolute-style: mean of per-CpG %-modified across covered CpGs in the window
mat_mean_pct <- build_matrix("mean_pct")
run_pca_plot(mat_mean_pct,
             "PCA of genome-wide 100kb methylation windows (mean %-modified)",
             "pca_windows.png")

# relative: methylated CpGs / covered CpGs per window, coverage-normalized
mat_relative <- build_matrix("relative_pct")
mat_relative_filtered <- mat_relative[apply(mat_relative, 1, var) > 0, ] # Identify and keep only the rows where variance is greater than 0
run_pca_plot(mat_relative_filtered,
             "PCA of genome-wide 100kb methylation windows (relative methylation)",
             "pca_windows_relative.png")

cat("Done. Plots written to:", normalizePath(plots_dir), "\n")
