# Metylace v genovych oblastech (gen + 3 kb upstream od zacatku genu) - PCA a rozdily vzorku
#
# Vstup: vystup genic_methylation_preprocess.sh (spustit nejdriv)
#   analyses/04_genic_methylation/output/<sample>_genes_up<UPSTREAM>_cov<MIN_COV>.tsv
#   chr start end gene_id gene_name strand gene_biotype n_valid n_mod n_cpg
#   metylace genu = sum(N_mod) / sum(N_valid) pres CpG 5mC volani v oblasti
#
# Vystup: matice gen x vzorek, PCA, korelace vzorku, rozdil pacientu vs. kontrola
#
# Pozadavky: R (data.table, ggplot2)
# Spusteni (pracovni adresar = methylation_profiles):
#   Rscript analyses/04_genic_methylation/genic_methylation_pca_diff.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- Parametry ----
SAMPLES    <- c("cfDNA_BAR", "cfDNA_MUD", "cfDNA_Pilot1", "cfDNA_VN")
CONTROL    <- "cfDNA_VN"                       # zdrava kontrola; ostatni = DLBCL

UPSTREAM   <- 3000          # musi odpovidat preprocess.sh (nazev vstupniho souboru)
MIN_COV    <- 5             # musi odpovidat preprocess.sh (nazev vstupniho souboru)
MIN_CPG    <- 10            # min. pocet CpG volani v oblasti (v kazdem vzorku)
TOP_VAR    <- 5000          # pocet nejvariabilnejsich genu do PCA
DELTA_MIN  <- 0.10          # min. |rozdil metylace| vuci kontrole

OUT_DIR  <- "analyses/04_genic_methylation/output"
PLOT_DIR <- "analyses/04_genic_methylation/plots"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Paleta / tema ----
CAT_COLORS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
SEQ_RAMP   <- c("#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b")
DIV_LOW <- "#184f95"; DIV_MID <- "#f0efec"; DIV_HIGH <- "#b23b3b"
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; SURFACE <- "#fcfcfb"

theme_meth <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(plot.background  = element_rect(fill = SURFACE, color = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = GRID, linewidth = 0.3),
          axis.text  = element_text(color = MUTED),
          axis.title = element_text(color = INK2),
          plot.title = element_text(color = INK, face = "bold"),
          plot.subtitle = element_text(color = INK2),
          legend.title  = element_text(color = INK2),
          legend.text   = element_text(color = INK2))
}
save_plot <- function(p, name, w, h) ggsave(file.path(PLOT_DIR, name), p, width = w, height = h, dpi = 200, bg = SURFACE)

# ---- Nacteni vystupu preprocess.sh ----
gene_counts_file <- function(s) file.path(OUT_DIR, sprintf("%s_genes_up%d_cov%d.tsv", s, UPSTREAM, MIN_COV))
missing <- SAMPLES[!file.exists(vapply(SAMPLES, gene_counts_file, ""))]
if (length(missing)) stop("Chybi vystup preprocess.sh pro: ", paste(missing, collapse = ", "),
                         "\nSpust nejdriv: bash analyses/04_genic_methylation/genic_methylation_preprocess.sh")

# ---- Matice gen x vzorek ----
read_counts <- function(s) {
  d <- fread(gene_counts_file(s), header = FALSE, na.strings = ".",
             col.names = c("chr", "start", "end", "gene_id", "gene_name", "strand", "gene_biotype",
                           "n_valid", "n_mod", "n_cpg"))
  d[, sample := s]
}
long <- rbindlist(lapply(SAMPLES, read_counts))
long[is.na(n_cpg), n_cpg := 0]
long[, meth := fifelse(n_valid > 0, n_mod / n_valid, NA_real_)]

global <- long[, .(n_genes_covered = sum(n_cpg >= MIN_CPG),
                   weighted_meth = sum(n_mod, na.rm = TRUE) / sum(n_valid, na.rm = TRUE),
                   median_gene_meth = median(meth[n_cpg >= MIN_CPG], na.rm = TRUE)), by = sample]
fwrite(global, file.path(OUT_DIR, "sample_summary.tsv"), sep = "\t")
print(global)

ok_genes <- long[, .(ok = all(n_cpg >= MIN_CPG & !is.na(meth)) && .N == length(SAMPLES)), by = gene_id][ok == TRUE, gene_id]
message("Genu s >= ", MIN_CPG, " CpG volanimi ve vsech vzorcich: ", length(ok_genes))

annot <- unique(long[, .(gene_id, gene_name, chr, start, end, strand)])
wide  <- dcast(long[gene_id %in% ok_genes], gene_id ~ sample, value.var = "meth")
wide  <- merge(annot, wide, by = "gene_id")
fwrite(wide, file.path(OUT_DIR, "genic_methylation_matrix.tsv"), sep = "\t")

mat <- as.matrix(wide[, ..SAMPLES]); rownames(mat) <- wide$gene_id

group <- setNames(ifelse(SAMPLES == CONTROL, "kontrola", "DLBCL"), SAMPLES)
GROUP_COLORS <- c(kontrola = CAT_COLORS[1], DLBCL = CAT_COLORS[2])

# ---- PCA ----
vars <- apply(mat, 1, var)
top  <- head(order(vars, decreasing = TRUE), min(TOP_VAR, nrow(mat)))
pca  <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)
ve   <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pc <- data.table(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                 PC3 = if (ncol(pca$x) >= 3) pca$x[, 3] else 0)
pc[, group := group[sample]]
fwrite(pc, file.path(OUT_DIR, "pca_scores.tsv"), sep = "\t")

load_dt <- data.table(gene_id = rownames(pca$rotation), PC1 = pca$rotation[, 1], PC2 = pca$rotation[, 2])
load_dt <- merge(annot[, .(gene_id, gene_name)], load_dt, by = "gene_id")[order(-abs(PC1))]
fwrite(load_dt, file.path(OUT_DIR, "pca_loadings.tsv"), sep = "\t")

pca_plot <- function(x, y, xi, yi) {
  ggplot(pc, aes(.data[[x]], .data[[y]], color = group)) +
    geom_hline(yintercept = 0, color = GRID) + geom_vline(xintercept = 0, color = GRID) +
    geom_point(size = 4) +
    geom_text(aes(label = sample), color = INK, vjust = -1.1, size = 3.8) +
    scale_color_manual(values = GROUP_COLORS, name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    labs(title = "PCA - metylace genovych oblasti",
         subtitle = sprintf("%d nejvariabilnejsich genu (gen + %d bp upstream), CpG 5mC, pokryti >= %d",
                            length(top), UPSTREAM, MIN_COV),
         x = sprintf("%s (%.1f %%)", x, ve[xi]), y = sprintf("%s (%.1f %%)", y, ve[yi])) +
    theme_meth() + theme(legend.position = "top", legend.justification = "left")
}
save_plot(pca_plot("PC1", "PC2", 1, 2), "pca_PC1_PC2.png", 7.5, 6.5)
if (ncol(pca$x) >= 3) save_plot(pca_plot("PC1", "PC3", 1, 3), "pca_PC1_PC3.png", 7.5, 6.5)

scree <- data.table(PC = factor(paste0("PC", seq_along(ve)), levels = paste0("PC", seq_along(ve))), ve = ve)
save_plot(ggplot(scree, aes(PC, ve)) + geom_col(fill = CAT_COLORS[1], width = 0.6) +
            geom_text(aes(label = sprintf("%.1f %%", ve)), vjust = -0.5, color = INK2) +
            labs(title = "Vysvetlena variabilita", x = NULL, y = "%") + theme_meth(),
          "pca_scree.png", 6, 4)

# ---- Korelace vzorku ----
cm <- cor(mat, method = "spearman")
cm_dt <- as.data.table(as.table(cm)); setnames(cm_dt, c("x", "y", "r"))
cm_dt[, `:=`(x = factor(x, levels = SAMPLES), y = factor(y, levels = rev(SAMPLES)))]
save_plot(ggplot(cm_dt, aes(x, y, fill = r)) +
            geom_tile(color = SURFACE, linewidth = 1.5) +
            geom_text(aes(label = sprintf("%.3f", r)), color = INK, size = 3.8) +
            scale_fill_gradientn(colors = SEQ_RAMP, name = "Spearman r") +
            labs(title = "Korelace vzorku - metylace genovych oblasti", x = NULL, y = NULL) +
            coord_fixed() + theme_meth() + theme(panel.grid = element_blank()),
          "sample_correlation.png", 7, 6)

# ---- Rozdily: pacienti vs. kontrola ----
patients <- setdiff(SAMPLES, CONTROL)
delta <- mat[, patients, drop = FALSE] - mat[, CONTROL]
diff_dt <- data.table(gene_id = rownames(mat),
                      control = mat[, CONTROL],
                      patient_mean = rowMeans(mat[, patients, drop = FALSE]),
                      delta_mean = rowMeans(delta),
                      n_hyper = rowSums(delta >= DELTA_MIN),
                      n_hypo  = rowSums(delta <= -DELTA_MIN))
diff_dt <- cbind(diff_dt, setnames(as.data.table(delta), paste0("delta_", patients)))
diff_dt[, direction := fcase(n_hyper == length(patients), "Hyper (vsichni pacienti)",
                             n_hypo  == length(patients), "Hypo (vsichni pacienti)",
                             default = "ostatni")]
diff_dt <- merge(annot, diff_dt, by = "gene_id")[order(-abs(delta_mean))]
fwrite(diff_dt, file.path(OUT_DIR, "genic_diff_vs_control.tsv"), sep = "\t")
fwrite(diff_dt[direction != "ostatni"], file.path(OUT_DIR, "genic_diff_vs_control_consistent.tsv"), sep = "\t")
message("Konzistentne hyper: ", sum(diff_dt$direction == "Hyper (vsichni pacienti)"),
        ", hypo: ", sum(diff_dt$direction == "Hypo (vsichni pacienti)"), " (|delta| >= ", DELTA_MIN, ")")

DIR_COLORS <- c("Hyper (vsichni pacienti)" = DIV_HIGH, "Hypo (vsichni pacienti)" = DIV_LOW, "ostatni" = "#c3c2b7")
lab_dt <- head(diff_dt[direction != "ostatni"], 15)

save_plot(ggplot(diff_dt[order(direction == "ostatni", decreasing = TRUE)],
                 aes(control, patient_mean, color = direction)) +
            geom_abline(slope = 1, intercept = c(-DELTA_MIN, 0, DELTA_MIN), color = MUTED,
                        linetype = c("dashed", "solid", "dashed"), linewidth = 0.3) +
            geom_point(size = 0.8, alpha = 0.6) +
            geom_text(data = lab_dt, aes(label = gene_name), color = INK, size = 3, vjust = -0.8, check_overlap = TRUE) +
            scale_color_manual(values = DIR_COLORS, name = NULL) +
            coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
            labs(title = "Metylace genovych oblasti: DLBCL vs. kontrola",
                 subtitle = sprintf("Kontrola = %s; barevne = |delta| >= %.2f ve vsech %d pacientech",
                                    CONTROL, DELTA_MIN, length(patients)),
                 x = paste0("Kontrola (", CONTROL, ")"), y = "Prumer DLBCL pacientu") +
            theme_meth() + theme(legend.position = "top", legend.justification = "left") +
            guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))),
          "diff_scatter_patients_vs_control.png", 7.5, 8)

delta_long <- melt(diff_dt[, c("gene_id", paste0("delta_", patients)), with = FALSE],
                   id.vars = "gene_id", variable.name = "sample", value.name = "delta")
delta_long[, sample := sub("^delta_", "", sample)]
save_plot(ggplot(delta_long, aes(delta, color = sample)) +
            geom_vline(xintercept = c(-DELTA_MIN, 0, DELTA_MIN), color = MUTED,
                       linetype = c("dashed", "solid", "dashed"), linewidth = 0.3) +
            geom_density(linewidth = 0.8) +
            scale_color_manual(values = setNames(CAT_COLORS[seq_along(patients) + 1], patients), name = NULL) +
            labs(title = paste0("Rozdil metylace vuci kontrole (", CONTROL, ")"),
                 subtitle = "Distribuce pres vsechny geny", x = "Delta metylace (pacient - kontrola)", y = "Hustota") +
            theme_meth() + theme(legend.position = "top", legend.justification = "left"),
          "diff_density.png", 8, 5)

# Heatmapa top rozdilnych genu (konzistentni smer, serazeno podle |delta|)
top_diff <- head(diff_dt[direction != "ostatni"], 50)
if (nrow(top_diff) == 0) top_diff <- head(diff_dt, 50)
hm <- melt(wide[gene_id %in% top_diff$gene_id, c("gene_id", "gene_name", SAMPLES), with = FALSE],
           id.vars = c("gene_id", "gene_name"), variable.name = "sample", value.name = "meth")
hm[, gene_label := factor(gene_name, levels = rev(unique(top_diff[order(delta_mean)]$gene_name)))]
hm[, sample := factor(sample, levels = c(CONTROL, patients))]
save_plot(ggplot(hm, aes(sample, gene_label, fill = meth)) +
            geom_tile(color = SURFACE, linewidth = 0.6) +
            scale_fill_gradientn(colors = c(DIV_MID, SEQ_RAMP), limits = c(0, 1),
                                 labels = scales::percent_format(accuracy = 1), name = "Metylace") +
            labs(title = "Geny s nejvetsim rozdilem DLBCL vs. kontrola",
                 subtitle = sprintf("Top %d podle |prumerne delta|; prvni sloupec = kontrola", nrow(top_diff)),
                 x = NULL, y = NULL) +
            theme_meth(10) +
            theme(panel.grid = element_blank(), axis.text.y = element_text(face = "italic", color = INK),
                  axis.text.x = element_text(color = INK, angle = 30, hjust = 1)),
          "diff_heatmap_top_genes.png", 6.5, max(5, 0.18 * nrow(top_diff) + 2))

message("Hotovo. Tabulky: ", OUT_DIR, "  Grafy: ", PLOT_DIR)
