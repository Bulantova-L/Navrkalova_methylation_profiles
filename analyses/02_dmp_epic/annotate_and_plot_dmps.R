# language: r
#
# Genome-wide plot of the differentially methylated EPIC probes found by
# differential_methylation_VN_vs_DLBCL.R, plus gene annotation of those sites.
#
# Input:  analyses/02_dmp_epic/output/dmp_significant.csv
#         (produced by differential_methylation_VN_vs_DLBCL.R)
#         data/hg38.chrom.sizes (for chromosome ordering/offsets)
# Run from methylation_profiles/: Rscript analyses/02_dmp_epic/annotate_and_plot_dmps.R
#
# Output: analyses/02_dmp_epic/plots/dmp_manhattan.png              - genome-wide plot of hits
#         analyses/02_dmp_epic/output/dmp_significant_annotated.csv
#                                               - hits + overlapping/nearest gene
#         analyses/02_dmp_epic/output/gene_hit_summary.csv
#                                               - genes ranked by number of hits
#
# Needs internet access (queries Ensembl via biomaRt), same approach already
# used elsewhere in this repo (extract_all_promoters.R).

library(data.table)
library(ggplot2)
library(biomaRt)
library(GenomicRanges)

dir.create("analyses/02_dmp_epic/output", recursive = TRUE, showWarnings = FALSE)
dir.create("analyses/02_dmp_epic/plots", recursive = TRUE, showWarnings = FALSE)

MAIN_CHROMS <- c(paste0("chr", 1:22), "chrX", "chrY")

# ---- Load DMP hits ----

dmp <- fread("analyses/02_dmp_epic/output/dmp_significant.csv")
dmp <- dmp[chr %in% MAIN_CHROMS]

# ---- Manhattan-style genome-wide plot ----

chrom_sizes <- fread("data/hg38.chrom.sizes", header = FALSE, col.names = c("chr", "size"))
chrom_sizes <- chrom_sizes[chr %in% MAIN_CHROMS]
chrom_sizes[, chr := factor(chr, levels = MAIN_CHROMS)]
setorder(chrom_sizes, chr)
chrom_sizes[, offset := cumsum(data.table::shift(size, fill = 0))]

dmp[, chr := factor(chr, levels = MAIN_CHROMS)]
dmp <- merge(dmp, chrom_sizes[, .(chr, offset)], by = "chr")
dmp[, genome_pos := start + offset]

chrom_mid <- chrom_sizes[, .(mid = offset + size / 2), by = chr]

manhattan_plot <- ggplot(dmp, aes(x = genome_pos, y = delta_beta, color = direction)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c(
    "Hyper (DLBCL > healthy)" = "#b23b3b",
    "Hypo (DLBCL < healthy)" = "#184f95"
  )) +
  scale_x_continuous(labels = chrom_mid$chr, breaks = chrom_mid$mid) +
  labs(
    x = "Chromosome", y = "Delta beta (mean DLBCL - healthy)", color = "Call",
    title = "Genome-wide distribution of differentially methylated EPIC probes"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))

ggsave("analyses/02_dmp_epic/plots/dmp_manhattan.png", plot = manhattan_plot, width = 14, height = 5, dpi = 300)

# ---- Gene annotation via Ensembl (GRCh38) ----

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 109)

genes <- as.data.table(getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "chromosome_name",
  values = c(1:22, "X", "Y"),
  mart = ensembl
))
genes <- genes[hgnc_symbol != ""]
genes[, chr := paste0("chr", chromosome_name)]

gene_gr <- GRanges(
  seqnames = genes$chr,
  ranges = IRanges(start = genes$start_position, end = genes$end_position),
  gene = genes$hgnc_symbol
)

dmp_gr <- GRanges(seqnames = dmp$chr, ranges = IRanges(start = dmp$start, end = dmp$end))

# Genes whose body overlaps the probe directly
overlap_hits <- findOverlaps(dmp_gr, gene_gr)
overlap_dt <- data.table(
  row = queryHits(overlap_hits),
  gene = gene_gr$gene[subjectHits(overlap_hits)]
)
overlap_genes <- overlap_dt[, .(overlapping_gene = paste(sort(unique(gene)), collapse = ";")), by = row]

# Nearest gene + distance, for every probe (useful when there is no direct overlap)
nearest_idx <- nearest(dmp_gr, gene_gr)
dmp[, nearest_gene := gene_gr$gene[nearest_idx]]
dmp[, nearest_gene_distance := distance(dmp_gr, gene_gr[nearest_idx])]

dmp[, row := .I]
dmp <- merge(dmp, overlap_genes, by = "row", all.x = TRUE)
dmp[, row := NULL]

fwrite(dmp, "analyses/02_dmp_epic/output/dmp_significant_annotated.csv")

# ---- Gene-level hit summary ----

gene_long <- dmp[!is.na(overlapping_gene), .(gene = unlist(strsplit(overlapping_gene, ";"))), by = .(probe_id, direction)]

gene_summary <- gene_long[, .(
  n_hits = .N,
  n_hyper = sum(direction == "Hyper (DLBCL > healthy)"),
  n_hypo = sum(direction == "Hypo (DLBCL < healthy)")
), by = gene]
setorder(gene_summary, -n_hits)

fwrite(gene_summary, "analyses/02_dmp_epic/output/gene_hit_summary.csv")

cat("Written:\n")
cat(" - analyses/02_dmp_epic/plots/dmp_manhattan.png\n")
cat(" - analyses/02_dmp_epic/output/dmp_significant_annotated.csv\n")
cat(" - analyses/02_dmp_epic/output/gene_hit_summary.csv\n")
cat("\nTop genes by number of differentially methylated probes:\n")
print(head(gene_summary, 20))
