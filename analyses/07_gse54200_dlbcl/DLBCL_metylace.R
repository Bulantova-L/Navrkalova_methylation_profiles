# Spoustet ze slozky analyses/07_gse54200_dlbcl (reference/ = GEO + liftOver, output/ = tabulky, plots/ = grafy)
setwd("/storage/713000-cmm/713016-bioit/base/workspace/lucka/Navrkalova_cfDNA/methylation_profiles/analyses/07_gse54200_dlbcl")
library(GEOquery)
gse <- getGEO("GSE54200")
class(gse)
gse1 <- if (is.list(gse)) gse[[1]] else gse

# metadata vzorků (skupiny DLBCL vs. NGCB) a anotace klonů/prob (souřadnice)
write.csv(pData(gse1), "reference/GSE54200_pheno.csv", row.names = TRUE)
write.csv(fData(gse1), "reference/GSE54200_features.csv", row.names = TRUE)

# samotná matice hodnot metylace (klony/proby × vzorky)
write.csv(exprs(gse1), "reference/GSE54200_matrix.csv", row.names = TRUE)


# Mann-Whitney test - určení odlišných míst -------------------------------

library(Biobase)
gse1 <- gse[[1]]
grp <- pData(gse1)$`cell type:ch1`
table(grp)

m <- exprs(gse1)  # sondy x vzorky
idx_dlbcl <- which(grp == "Diffuse Large B cell Lymphoma cells")
idx_cb    <- which(grp == "normal Germinal Center B Cells")

pvals <- apply(m, 1, function(x) {
  a <- x[idx_dlbcl]; b <- x[idx_cb]
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) >= 5 && length(b) >= 5) wilcox.test(a, b)$p.value else NA
})

med_dlbcl <- apply(m[, idx_dlbcl], 1, median, na.rm = TRUE)
med_cb    <- apply(m[, idx_cb], 1, median, na.rm = TRUE)

feat <- fData(gse1)
result <- data.frame(
  probe_id = rownames(m),
  chrom = feat$CHROMOSOME, start = feat$RANGE_START, end = feat$RANGE_STOP,
  DLBCL_median = med_dlbcl, CB_median = med_cb, delta = med_dlbcl - med_cb,
  p_value = pvals, fdr = p.adjust(pvals, method = "BH")
)
result <- result[order(result$p_value), ]
write.csv(result, "output/GSE54200_DLBCLvsCB_mannwhitney.csv", row.names = FALSE)

### LIFT OVER to hg38
library(rtracklayer)
# chain soubor pro import.chain musí být rozbalený (ne .gz)
if (!file.exists("reference/hg17ToHg38.over.chain")) system("gunzip -k reference/hg17ToHg38.over.chain.gz")
chain <- import.chain("reference/hg17ToHg38.over.chain")

gr17 <- GRanges(seqnames = result$chrom,
                ranges = IRanges(start = result$start, end = result$end),
                probe_id = result$probe_id)
gr38 <- unlist(liftOver(gr17, chain))

# export do BED (0-based, half-open) pro bedtools intersect
gr38_bed <- data.frame(
  chrom = as.character(seqnames(gr38)),
  start = start(gr38) - 1L,
  end   = end(gr38),
  name  = gr38$probe_id
)
write.table(gr38_bed, "reference/output_hg38_new.bed", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)


# Analysis on our samples -------------------------------------------------

library(dplyr)
library(purrr)

cols <- c("chrom","start","end","probe_id","Nmod","Ncov","pct_modified")

sample_ids <- c("BAR", "MUD", "Pilot1", "VN")
samples <- setNames(
  file.path("..", "..", "results", paste0("cfDNA_", sample_ids),   # Snakefile rule bedtools_gse54200_map
            "methylation", paste0("cfDNA_", sample_ids, "_gse54200_mapped.bed")),
  paste0("cfDNA_", sample_ids)
)

control_sample <- "cfDNA_VN"   
patient_samples <- setdiff(names(samples), control_sample)

agg_list <- lapply(samples, function(path) {
  df <- read.table(path, sep = "\t", header = FALSE, col.names = cols,
                   na.strings = "NA", stringsAsFactors = FALSE)
  df %>%
    group_by(probe_id) %>%
    summarise(Nmod = sum(Nmod, na.rm = TRUE),
              Ncov = sum(Ncov, na.rm = TRUE), .groups = "drop")
})

names(agg_list) <- names(samples)

wide <- Reduce(function(acc, nm) {
  d <- agg_list[[nm]] %>% rename_with(~paste0(.x, "_", nm), c(Nmod, Ncov))
  if (is.null(acc)) d else full_join(acc, d, by = "probe_id")
}, names(agg_list), NULL)

lit <- read.csv("output/GSE54200_DLBCLvsCB_mannwhitney.csv") %>%
  dplyr::rename(delta_DLBCL_minus_CB = delta, fdr_bh = fdr)  # tvoje vlastní R reprodukce má jiné názvy sloupců než můj původní soubor

fisher_per_patient <- function(patient) {
  npt <- wide[[paste0("Nmod_", patient)]]; cpt <- wide[[paste0("Ncov_", patient)]]
  nc  <- wide[[paste0("Nmod_", control_sample)]]; cc <- wide[[paste0("Ncov_", control_sample)]]
  pmap_dfr(list(wide$probe_id, npt, cpt, nc, cc), function(probe_id, npt, cpt, nc, cc) {
    if (any(is.na(c(npt, cpt, nc, cc))) || cpt == 0 || cc == 0) return(NULL)
    tab <- matrix(c(npt, cpt - npt, nc, cc - nc), nrow = 2)
    p <- fisher.test(tab)$p.value
    tibble::tibble(patient = patient, probe_id = probe_id,
                   pct_pt = 100 * npt / cpt, pct_control = 100 * nc / cc,
                   delta = 100 * npt / cpt - 100 * nc / cc, p_value = p)
  })
}

res <- bind_rows(lapply(patient_samples, fisher_per_patient)) %>%
  group_by(patient) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  left_join(lit %>% select(probe_id, delta_DLBCL_minus_CB, fdr_bh_literatura = fdr_bh),
            by = "probe_id")

write.csv(res, "output/vysledky_pacienti_vs_kontrola.csv", row.names = FALSE)

sig <- res %>% filter(fdr < 0.05) %>%
  mutate(smer_shoduje = (delta * delta_DLBCL_minus_CB * -1) > 0)
sig %>% group_by(patient) %>% summarise(podil_shody = mean(smer_shoduje, na.rm = TRUE))


sig %>%
  group_by(patient) %>%
  summarise(
    n_loci = n(),
    n_shoda = sum(smer_shoduje, na.rm = TRUE),
    podil_shody = mean(smer_shoduje, na.rm = TRUE),
    p_binom = binom.test(n_shoda, n_loci, p = 0.5, alternative = "greater")$p.value
  )

# probe_id nemusi byt po liftOveru jedinecny klic: jeden hg17 lokus se pri liftOveru
# muze rozpadnout do vice hg38 useku, ktere si podrzi stejne probe_id ale ruzne
# souradnice/hodnoty -> zkontroluj, jestli k tomu u tebe dochazi:
res %>% count(patient, probe_id) %>% filter(n > 1) %>% nrow()

# pro kazdou dvojici (patient, probe_id) nechej jen jeden radek - ten s nejnizsi p-hodnotou
# (nejspolehlivejsi odhad), jinak joiny padaji na duplicitach
res <- res %>%
  group_by(patient, probe_id) %>%
  slice_min(p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

# loci signifikantni aspon u jednoho pacienta (FDR < 0.05) - jinak je heatmapa neprehledna
sig_probes <- res %>% filter(fdr < 0.05) %>% distinct(probe_id) %>% pull(probe_id)
length(sig_probes)  # kolik jich celkem je - pokud stovky+, heatmapa dole se omezi na TOP_N nejsilnejsich

# kontrola (pct_control, delta_DLBCL_minus_CB) je v ramci jednoho probe_id stejna,
# ale kdyby i literaturni tabulka mela na probe_id vic radku, znovu vezmi jen jeden
control_vals <- res %>%
  filter(probe_id %in% sig_probes) %>%
  distinct(probe_id, pct_control, delta_DLBCL_minus_CB) %>%
  group_by(probe_id) %>%
  slice_min(abs(delta_DLBCL_minus_CB), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(smer_literatura = ifelse(delta_DLBCL_minus_CB < 0,
                                  "hypermetylace v DLBCL", "hypometylace v DLBCL")) %>%
  arrange(smer_literatura, desc(abs(delta_DLBCL_minus_CB)))

# heatmapa je citelna jen do ~100-200 radku (a png/cairo ma i tvrdy strop na
# rozmer obrazku) - omez zobrazeni na TOP_N nejsilnejsich lokusu v kazdem smeru
# (hyper/hypo) zvlast. Kompletni tabulka vsech signifikantnich lokusu zustava v
# `control_vals` / `vysledky_pacienti_vs_kontrola.csv`, jen se nezobrazi cela v obrazku.
TOP_N_PER_SMER <- 75
control_vals_plot <- control_vals %>%
  group_by(smer_literatura) %>%
  slice_head(n = TOP_N_PER_SMER) %>%
  ungroup()

# HEATMAP -----------------------------------------------------------------

library(pheatmap)
library(tidyr)
library(tibble)
res <- read.csv("output/vysledky_pacienti_vs_kontrola.csv")

mat_patients <- res %>%
  dplyr::filter(probe_id %in% control_vals_plot$probe_id) %>%
  dplyr::select(patient, probe_id, pct_pt) %>%
  tidyr::pivot_wider(names_from = patient, values_from = pct_pt,
                     values_fn = ~ mean(.x, na.rm = TRUE)) %>%  # pojistka, kdyby v res porad byly duplicity
  tibble::column_to_rownames("probe_id")

mat <- as.matrix(mat_patients)
mat <- cbind(mat, "cfDNA_VN (kontrola)" = control_vals_plot$pct_control[match(rownames(mat), control_vals_plot$probe_id)])
mat <- mat[control_vals_plot$probe_id, , drop = FALSE]  # poradi radku podle smeru z literatury
storage.mode(mat) <- "numeric"  # dalsi pojistka - kdyby cbind/pivot_wider necim zmenilo typ na list/character

row_anno <- data.frame(smer_literatura = control_vals_plot$smer_literatura,
                       row.names = control_vals_plot$probe_id)

pheatmap(mat,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = row_anno,
         color = viridisLite::viridis(100),
         breaks = seq(0, 100, length.out = 101),
         fontsize_row = 6,
         main = paste0("Metylace na top ", nrow(mat), " signifikantnich loci (z ", length(sig_probes),
                       " celkem, FDR < 0.05) vs. smer z GSE54200"),
         filename = "plots/heatmapa_metylace.png", width = 8,
         height = min(30, max(4, 0.15 * nrow(mat))))  # strop 30" - png/cairo ma tvrdy limit na rozmer obrazku
