# language: r
#
# Differentially methylated positions (DMPs): healthy control (cfDNA_VN)
# vs. DLBCL patients (cfDNA_BAR, cfDNA_LOS, cfDNA_MUD, cfDNA_Pilot1).
#
# Input:  results/<sample>/methylation/<sample>_epic_mapped.bed
#         (chr, start, end, probe_id, methylation % [0-100] per EPIC probe site)
# Method: delta-beta threshold on mean(patients) - control, per probe.
#         (n = 1 control sample -> no group variance available, so this is a
#         magnitude-based call, not a p-value-based one.)
# Run from methylation_profiles/: Rscript analyses/02_dmp_epic/differential_methylation_VN_vs_DLBCL.R

library(data.table)
library(ggplot2)

CONTROL_SAMPLE <- "cfDNA_VN"
DELTA_THRESHOLD <- 20   # percentage points (i.e. beta-value difference of 0.2)

dir.create("analyses/02_dmp_epic/output", recursive = TRUE, showWarnings = FALSE)
dir.create("analyses/02_dmp_epic/plots", recursive = TRUE, showWarnings = FALSE)

# ---- Load ----

read_epic_file <- function(file_path) {
  sample_id <- sub("_epic_mapped\\.bed$", "", basename(file_path))

  dt <- fread(file_path, header = FALSE, na.strings = ".")
  setnames(dt, c("chr", "start", "end", "probe_id", "methylation"))
  dt[, sample := sample_id]
  dt[, methylation := as.numeric(methylation)]
  dt
}

files <- list.files("results", pattern = "_epic_mapped\\.bed$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No *_epic_mapped.bed files found under results/.")

all_data <- rbindlist(lapply(files, read_epic_file))

samples <- sort(unique(all_data$sample))
patient_samples <- setdiff(samples, CONTROL_SAMPLE)
if (!(CONTROL_SAMPLE %in% samples)) stop("Control sample '", CONTROL_SAMPLE, "' not found among: ", paste(samples, collapse = ", "))

cat("Samples found:", paste(samples, collapse = ", "), "\n")
cat("Control:", CONTROL_SAMPLE, "| Patients:", paste(patient_samples, collapse = ", "), "\n")

# ---- Wide matrix (probe x sample) ----

meth_wide <- dcast(all_data, chr + start + end + probe_id ~ sample, value.var = "methylation")

n_total <- nrow(meth_wide)

# Complete-case filter: require methylation at the control and every patient sample
complete_rows <- complete.cases(meth_wide[, c(CONTROL_SAMPLE, patient_samples), with = FALSE])
meth_complete <- meth_wide[complete_rows]

cat("Probes total:", n_total, "\n")
cat("Probes with complete data across all samples:", nrow(meth_complete), "\n")
cat("Probes dropped (missing coverage in >=1 sample):", n_total - nrow(meth_complete), "\n")

# ---- Delta-beta ----

meth_complete[, control_meth := get(CONTROL_SAMPLE)]
meth_complete[, patient_mean_meth := rowMeans(.SD), .SDcols = patient_samples]
meth_complete[, delta_beta := patient_mean_meth - control_meth]

meth_complete[, direction := fcase(
  delta_beta >= DELTA_THRESHOLD, "Hyper (DLBCL > healthy)",
  delta_beta <= -DELTA_THRESHOLD, "Hypo (DLBCL < healthy)",
  default = "NS"
)]

meth_complete <- meth_complete[order(-abs(delta_beta))]

# ---- Outputs ----

fwrite(meth_complete, "analyses/02_dmp_epic/output/dmp_results.csv")

summary_counts <- meth_complete[, .N, by = direction]
cat("\nClassification summary (|delta beta| >=", DELTA_THRESHOLD, "percentage points):\n")
print(summary_counts)

top_hits <- meth_complete[direction != "NS"]
fwrite(top_hits, "analyses/02_dmp_epic/output/dmp_significant.csv")

n_hyper <- summary_counts[direction == "Hyper (DLBCL > healthy)", N]
n_hypo <- summary_counts[direction == "Hypo (DLBCL < healthy)", N]
cat("\nHypermethylated in DLBCL:", if (length(n_hyper)) n_hyper else 0, "\n")
cat("Hypomethylated in DLBCL:", if (length(n_hypo)) n_hypo else 0, "\n")

# ---- Plots ----

dmp_scatter <- ggplot(meth_complete, aes(x = control_meth, y = patient_mean_meth, color = direction)) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_abline(slope = 1, intercept = DELTA_THRESHOLD, linetype = "dashed", color = "grey40") +
  geom_abline(slope = 1, intercept = -DELTA_THRESHOLD, linetype = "dashed", color = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
  scale_color_manual(values = c(
    "Hyper (DLBCL > healthy)" = "#b23b3b",
    "Hypo (DLBCL < healthy)" = "#184f95",
    "NS" = "grey80"
  )) +
  labs(
    x = paste0(CONTROL_SAMPLE, " methylation (%)"),
    y = "Mean DLBCL patient methylation (%)",
    color = "Call",
    title = "EPIC-probe methylation: healthy control vs. DLBCL patients"
  ) +
  theme_minimal(base_size = 13)

ggsave("analyses/02_dmp_epic/plots/dmp_scatter.png", plot = dmp_scatter, width = 8, height = 7, dpi = 300)

dmp_hist <- ggplot(meth_complete, aes(x = delta_beta)) +
  geom_histogram(bins = 100, fill = "#4a3aa7") +
  geom_vline(xintercept = c(-DELTA_THRESHOLD, DELTA_THRESHOLD), linetype = "dashed", color = "grey40") +
  labs(
    x = "Delta beta (mean DLBCL - healthy), percentage points",
    y = "Number of probes",
    title = "Distribution of methylation differences per EPIC probe"
  ) +
  theme_minimal(base_size = 13)

ggsave("analyses/02_dmp_epic/plots/dmp_histogram.png", plot = dmp_hist, width = 8, height = 5, dpi = 300)

cat("\nWritten:\n")
cat(" - analyses/02_dmp_epic/output/dmp_results.csv (all tested probes)\n")
cat(" - analyses/02_dmp_epic/output/dmp_significant.csv (hyper/hypo only)\n")
cat(" - analyses/02_dmp_epic/plots/dmp_scatter.png\n")
cat(" - analyses/02_dmp_epic/plots/dmp_histogram.png\n")
