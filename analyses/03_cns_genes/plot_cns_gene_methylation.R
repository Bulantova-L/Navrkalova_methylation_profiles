# Vizualizace metylace CNS genu (promotory, TSS -2kb/+500bp) napric vzorky
#
# Vstup:  results/<sample>/methylation/<sample>_cns_mapped.bed
#         (chr, start, end, gene, prumer sloupce 11 modkit = % modifikace, pocet zaznamu) - viz rule bedtools_cns_aggregate
#
# POZOR: modkit pileup bezel bez --cpg, takze bedtools map prumeruje % modifikace pres VSECHNY C
# (CpG i non-CpG, oba retezce, kody m i h). Hodnoty jsou proto v procentech a "zredene" (0-~13 %),
# nejde o klasickou CpG beta hodnotu. Pro relativni srovnani vzorku/genu je to pouzitelne.
# Vystup: analyses/03_cns_genes/plots/cns_gene_methylation_heatmap.png   - gen x vzorek, vazeny prumer
#         analyses/03_cns_genes/plots/cns_gene_methylation_profiles.png  - profil podel genomu pro kazdy gen
#         analyses/03_cns_genes/output/cns_gene_methylation_summary.tsv
#
# Spusteni (pracovni adresar = methylation_profiles):
#   Rscript analyses/03_cns_genes/plot_cns_gene_methylation.R

library(data.table)
library(ggplot2)

# ---- Paleta (reference palette, light mode) ----
CAT_COLORS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
SEQ_RAMP   <- c("#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b")
GRID_COLOR <- "#e1e0d9"
AXIS_COLOR <- "#c3c2b7"
MUTED_COLOR <- "#898781"
SECONDARY_INK <- "#52514e"
INK_COLOR  <- "#0b0b0b"
SURFACE    <- "#fcfcfb"

theme_meth <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background  = element_rect(fill = SURFACE, color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = GRID_COLOR, linewidth = 0.3),
      axis.text        = element_text(color = MUTED_COLOR),
      axis.title       = element_text(color = SECONDARY_INK),
      plot.title       = element_text(color = INK_COLOR, face = "bold"),
      plot.subtitle    = element_text(color = SECONDARY_INK),
      legend.title     = element_text(color = SECONDARY_INK),
      legend.text      = element_text(color = SECONDARY_INK),
      strip.text       = element_text(color = INK_COLOR, face = "bold", hjust = 0)
    )
}

# ---- Nacteni dat ----
bed_files <- list.files("results", pattern = "_cns_mapped\\.bed$", recursive = TRUE, full.names = TRUE)
if (length(bed_files) == 0) stop("Nenalezen zadny results/*/methylation/*_cns_mapped.bed")

read_cns <- function(f) {
  dt <- fread(f, header = FALSE, col.names = c("chr", "start", "end", "gene", "methylation", "n_calls"),
              na.strings = c(".", "NA"))
  dt[, sample := sub("_cns_mapped\\.bed$", "", basename(f))]
  unique(dt)  # promotorovy BED obsahuje duplicitni radky (vice transkriptu se stejnym TSS)
}
d <- rbindlist(lapply(bed_files, read_cns))
d[, methylation := as.numeric(methylation) / 100]  # % -> podil
d[, n_calls := as.numeric(n_calls)]
d[, pos_mb := (start + end) / 2 / 1e6]

samples <- sort(unique(d$sample))
sample_colors <- setNames(CAT_COLORS[seq_along(samples)], samples)
sample_shapes <- setNames(c(16, 17, 15, 18, 3, 4, 8, 1)[seq_along(samples)], samples)

# ---- Souhrn gen x vzorek (vazeny poctem CpG volani) ----
gene_sum <- d[!is.na(methylation) & n_calls > 0,
              .(methylation = sum(methylation * n_calls) / sum(n_calls),  # vazeno poctem zaznamu
                n_windows = .N, n_calls = sum(n_calls)),
              by = .(gene, sample)]
dir.create("analyses/03_cns_genes/output", recursive = TRUE, showWarnings = FALSE)
fwrite(dcast(gene_sum, gene ~ sample, value.var = "methylation"),
       "analyses/03_cns_genes/output/cns_gene_methylation_summary.tsv", sep = "\t")

# Geny serazene podle prumerne metylace napric vzorky
gene_order <- gene_sum[, .(m = mean(methylation)), by = gene][order(m), gene]
gene_sum[, gene := factor(gene, levels = gene_order)]
gene_sum[, sample := factor(sample, levels = samples)]

heat_max <- ceiling(max(gene_sum$methylation) * 100) / 100

p_heat <- ggplot(gene_sum, aes(x = sample, y = gene, fill = methylation)) +
  geom_tile(color = SURFACE, linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.1f", 100 * methylation),
                color = methylation > 0.55 * heat_max), size = 3.6, show.legend = FALSE) +
  scale_color_manual(values = c(`FALSE` = INK_COLOR, `TRUE` = "#ffffff")) +
  scale_fill_gradientn(colors = SEQ_RAMP, limits = c(0, heat_max),
                       labels = scales::percent_format(accuracy = 1), name = "Modifikace") +
  labs(title = "Metylace promotoru CNS genu v cfDNA",
       subtitle = "Prumerne % modifikace (5mC+5hmC, vsechna C) v promotorech TSS -2 kb/+500 bp\nHodnota v bunce = %; vazeno poctem zaznamu; prazdne = bez pokryti",
       x = NULL, y = NULL) +
  theme_meth() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = INK_COLOR, face = "italic"),
        axis.text.x = element_text(color = INK_COLOR),
        legend.key.height = unit(1.6, "cm"))

ggsave("analyses/03_cns_genes/plots/cns_gene_methylation_heatmap.png", p_heat, width = 7.5, height = 7, dpi = 200, bg = SURFACE)

# ---- Profil podel genomu pro kazdy gen ----
d[, gene := factor(gene, levels = rev(gene_order))]
d[, sample := factor(sample, levels = samples)]
d[, ref_pos := min(start + end) / 2, by = gene]
d[, pos_kb := ((start + end) / 2 - ref_pos) / 1e3]
d[, panel := sprintf("%s  (%s:%s)", gene, chr, format(ref_pos, big.mark = ",", scientific = FALSE, trim = TRUE))]
panel_levels <- unique(d[order(gene), panel])
d[, panel := factor(panel, levels = panel_levels)]

p_prof <- ggplot(d[!is.na(methylation)],
                 aes(x = pos_kb, y = methylation, color = sample, shape = sample, group = sample)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 2.4) +
  facet_wrap(~ panel, scales = "free_x", ncol = 4) +
  scale_color_manual(values = sample_colors, name = "Vzorek") +
  scale_shape_manual(values = sample_shapes, name = "Vzorek") +
  scale_y_continuous(limits = c(0, NA), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.08, add = 0.5)) +
  labs(title = "Metylacni profil promotoru CNS genu",
       subtitle = "Kazdy bod = jedno promotorove okno (alternativni TSS); geny serazeny od nejvice metylovaneho",
       x = "Vzdalenost od prvniho promotoroveho okna genu (kb, hg38)", y = "Modifikace (% pres vsechna C)") +
  theme_meth(11) +
  theme(legend.position = "top", legend.justification = "left",
        panel.border = element_rect(color = AXIS_COLOR, fill = NA, linewidth = 0.3))

ggsave("analyses/03_cns_genes/plots/cns_gene_methylation_profiles.png", p_prof, width = 13, height = 10, dpi = 200, bg = SURFACE)

cat("Hotovo:", nrow(gene_sum), "gen x vzorek kombinaci,", nrow(d), "oken\n")
