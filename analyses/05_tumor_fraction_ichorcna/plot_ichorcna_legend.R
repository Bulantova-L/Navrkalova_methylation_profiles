# Genome-wide CNV grafy z vystupu ichorCNA - s legendou stavu kopiovosti
#
# Vstup:  analyses/05_tumor_fraction_ichorcna/output/ichorcna*/<pon>/<sample>/<sample>.cna.seg   (biny: logR + Corrected_Call)
#         analyses/05_tumor_fraction_ichorcna/output/ichorcna*/<pon>/<sample>/<sample>.seg.txt   (segmenty: median logR + Corrected_Call)
#         analyses/05_tumor_fraction_ichorcna/output/ichorcna*/<pon>/<sample>/<sample>.params.txt (tumor fraction, ploidie)
# Vystup: analyses/05_tumor_fraction_ichorcna/output/ichorcna*/plots_legend/<pon>_<sample>_genomeWide.png  (jeden vzorek)
#         analyses/05_tumor_fraction_ichorcna/output/ichorcna*/plots_legend/<pon>_all_samples_genomeWide.png (vsechny vzorky pod sebou)
#
# Barvy odpovidaji puvodnim grafum ichorCNA (USE_ICHOR_COLORS = TRUE),
# nebo divergentni skala ztrata = modra / zisk = cervena (FALSE).
#
# Spusteni (pracovni adresar = methylation_profiles):
#   Rscript analyses/05_tumor_fraction_ichorcna/plot_ichorcna_legend.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- Parametry ----
ICHOR_DIRS <- Sys.glob("analyses/05_tumor_fraction_ichorcna/output/ichorcna*")    # vsechny behy (plne delky i frag90-150 apod.)
USE_ICHOR_COLORS <- TRUE
OUT_OVERRIDE <- Sys.getenv("ICHOR_PLOT_DIR", "")   # prazdne = <beh>/plots_legend; jinak spolecna slozka
CHRS <- paste0("chr", c(1:22, "X"))

CALL_LEVELS <- c("HOMD", "HETD", "NEUT", "GAIN", "AMP", "HLAMP")
CALL_LABELS <- c(HOMD  = "HOMD - homozygous deletion (0 copies)",
                 HETD  = "HETD - loss (1 copy)",
                 NEUT  = "NEUT - neutral (2 copies)",
                 GAIN  = "GAIN - gain (3 copies)",
                 AMP   = "AMP - amplification (4 copies)",
                 HLAMP = "HLAMP - high-level amplification (5+ copies)")
CALL_COLORS <- if (USE_ICHOR_COLORS) {
  c(HOMD = "#00FF00", HETD = "#006400", NEUT = "#0000FF", GAIN = "#8B0000", AMP = "#FF0000", HLAMP = "#FF0000")
} else {
  c(HOMD = "#0d366b", HETD = "#3987e5", NEUT = "#9e9d97", GAIN = "#e66767", AMP = "#b23b3b", HLAMP = "#7a1f1f")
}

INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; SURFACE <- "#fcfcfb"; BAND <- "#f0efec"

# ---- Pomocne funkce ----
norm_call <- function(x) {
  x <- sub("^HLAMP[0-9]+$", "HLAMP", x)
  factor(x, levels = CALL_LEVELS)
}

read_params <- function(f) {
  p <- readLines(f, n = 2)
  v <- strsplit(p[2], "\t")[[1]]
  list(tf = as.numeric(v[2]), ploidy = as.numeric(v[3]))
}

load_sample <- function(dir) {
  s <- basename(dir)
  cna_f <- file.path(dir, paste0(s, ".cna.seg"))
  seg_f <- file.path(dir, paste0(s, ".seg.txt"))
  par_f <- file.path(dir, paste0(s, ".params.txt"))
  if (!all(file.exists(c(cna_f, seg_f, par_f)))) return(NULL)

  cna <- fread(cna_f)
  setnames(cna, sub(paste0("^", s, "\\."), "", names(cna)))   # "cfDNA_X.logR" -> "logR"
  cna <- cna[chr %in% CHRS & !is.na(logR), .(chr, start, end, logR, call = norm_call(Corrected_Call))]

  seg <- fread(seg_f)[chrom %in% CHRS, .(chr = chrom, start, end, logR = seg.median.logR,
                                          call = norm_call(Corrected_Call))]
  prm <- read_params(par_f)
  list(sample = s, cna = cna, seg = seg, tf = prm$tf, ploidy = prm$ploidy)
}

# Kumulativni genomove souradnice
chrom_layout <- function(sizes_file) {
  cs <- fread(sizes_file, header = FALSE, col.names = c("chr", "len"))[chr %in% CHRS]
  cs <- cs[match(CHRS, chr)][!is.na(chr)]
  cs[, offset := c(0, cumsum(as.numeric(len))[-.N])]
  cs[, mid := offset + len / 2]
  cs[, label := sub("^chr", "", chr)]
  cs[, band := seq_len(.N) %% 2 == 0]
  cs
}

to_genome <- function(dt, lay) {
  dt <- merge(dt, lay[, .(chr, offset)], by = "chr")
  dt[, `:=`(gstart = start + offset, gend = end + offset, gmid = (start + end) / 2 + offset)]
  dt
}

plot_samples <- function(samples, lay, pon) {
  cna <- rbindlist(lapply(samples, function(x) to_genome(x$cna, lay)[, sample := x$sample]))
  seg <- rbindlist(lapply(samples, function(x) to_genome(x$seg, lay)[, sample := x$sample]))
  lab <- vapply(samples, function(x) sprintf("%s   (TF %.1f %%, ploidy %.2f)", x$sample, 100 * x$tf, x$ploidy), "")
  names(lab) <- vapply(samples, `[[`, "", "sample")
  cna[, panel := factor(lab[sample], levels = lab)]
  seg[, panel := factor(lab[sample], levels = lab)]

  ylim <- c(-2, 2)
  cna[, logR_plot := pmin(pmax(logR, ylim[1]), ylim[2])]   # hodnoty mimo rozsah na okraj (jako ichorCNA)

  ggplot() +
    geom_rect(data = lay[band == TRUE], aes(xmin = offset, xmax = offset + len, ymin = -Inf, ymax = Inf),
              fill = BAND, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, color = MUTED, linewidth = 0.3) +
    geom_point(data = cna, aes(gmid, logR_plot, color = call), size = 0.6, alpha = 0.8) +
    geom_segment(data = seg, aes(x = gstart, xend = gend, y = logR, yend = logR, color = call), linewidth = 1.1) +
    facet_wrap(~ panel, ncol = 1) +
    scale_color_manual(values = CALL_COLORS, labels = CALL_LABELS, breaks = CALL_LEVELS,
                       drop = FALSE, name = "Copy-number state (Corrected_Call)") +
    scale_x_continuous(breaks = lay$mid, labels = lay$label, expand = expansion(mult = 0.005)) +
    scale_y_continuous(limits = ylim, breaks = -2:2) +
    labs(title = paste0("ichorCNA - genome-wide copy number (", pon, ")"),
         subtitle = "Points = 1 Mb bins, lines = segments; log2 ratio vs. panel of healthy normals",
         x = "Chromosome", y = "Copy number (log2 ratio)") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1, linewidth = 0), nrow = 2)) +
    theme_minimal(base_size = 11) +
    theme(plot.background  = element_rect(fill = SURFACE, color = NA),
          panel.grid       = element_blank(),
          panel.grid.major.y = element_line(color = GRID, linewidth = 0.3),
          axis.text        = element_text(color = MUTED),
          axis.title       = element_text(color = INK2),
          strip.text       = element_text(color = INK, face = "bold", hjust = 0),
          plot.title       = element_text(color = INK, face = "bold"),
          plot.subtitle    = element_text(color = INK2),
          legend.position  = "top", legend.justification = "left",
          legend.title     = element_text(color = INK2),
          legend.text      = element_text(color = INK2))
}

# ---- Hlavni smycka ----
for (idir in ICHOR_DIRS) {
  sizes <- file.path(idir, "chrom_sizes.txt")
  if (!file.exists(sizes)) next
  lay <- chrom_layout(sizes)
  out_dir <- if (nzchar(OUT_OVERRIDE)) file.path(OUT_OVERRIDE, basename(idir)) else file.path(idir, "plots_legend")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Nelze vytvorit ", out_dir, " (prava?) - nastav ICHOR_PLOT_DIR")

  for (pdir in Sys.glob(file.path(idir, "pon_*"))) {
    pon <- basename(pdir)
    samples <- Filter(Negate(is.null), lapply(list.dirs(pdir, recursive = FALSE), load_sample))
    if (!length(samples)) next

    for (x in samples) {
      ggsave(file.path(out_dir, sprintf("%s_%s_genomeWide.png", pon, x$sample)),
             plot_samples(list(x), lay, pon), width = 14, height = 4.2, dpi = 200, bg = SURFACE)
    }
    ggsave(file.path(out_dir, sprintf("%s_all_samples_genomeWide.png", pon)),
           plot_samples(samples, lay, pon), width = 14, height = 1.2 + 2.3 * length(samples),
           dpi = 200, bg = SURFACE, limitsize = FALSE)
    message(idir, " / ", pon, ": ", length(samples), " vzorku -> ", out_dir)
  }
}
