# language: r
# Potřebné balíčky: tidyverse, limma
library(tidyverse)
library(limma)

# 1) Sestav probe x sample z vašich matched_cpgs.bed souborů (pokud máte 5 souborů)
files <- dir("methylation_profiles", pattern = "matched_cpgs.*\\.bed$", full.names = TRUE)
# pokud máte jinou konvenci jmen, upravte pattern nebo použijte explicitní vektor cest

sample_tbls <- map(files, ~ read_tsv(.x, col_names = FALSE, show_col_types = FALSE) %>%
  transmute(probe = X22,      # vaše X22 obsahuje "cgID_suffix"
            beta = X11 / 100) %>%   # X11 percent -> 0-1 beta
  mutate(probe = sub("_.*$", "", probe)) %>% # odstranit _BC11 apod -> cgID
  group_by(probe) %>% summarize(beta = mean(beta, na.rm = TRUE)) )

names(sample_tbls) <- basename(files) %>% tools::file_path_sans_ext()

beta_mat <- reduce(sample_tbls, full_join, by = "probe") %>%
  arrange(probe) %>%
  replace_na(list()) # ensure no missing tibble structure
# přejmenuj sloupce a vytvoř matici
colnames(beta_mat)[-1] <- names(sample_tbls)
beta_m <- as.matrix(beta_mat[,-1])
rownames(beta_m) <- beta_mat$probe

# 2) Načíst referenční EPIC matrix (řádky = cgIDs, sloupce = referenční vzorky), a labely
# předpoklad: reference_epic_matrix.txt má řádky = probe IDs (rownames), sloupce = vzorky, hodnoty = beta(0-1)
ref <- read_tsv("reference_epic_matrix.txt", col_names = TRUE) # upravte pokud jiný formát
# pokud soubor má rownames v první kolonce:
# ref <- read_tsv("reference_epic_matrix.txt", col_names = TRUE) %>% column_to_rownames(var = "probe")

# Zjednodušeně: vytvořte data.frame ref_beta (rows probes, cols samples) a ref_labels vektor s typy (PCNSL/systemic/normal)
# --- UŽIVATEL musí zajistit ref_beta a ref_labels ---
# Příklad (nahraďte skutečným načtením):
# ref_beta <- as.matrix(ref)
# ref_labels <- c("PCNSL","other","normal", ...)

# 3) Propojit průnik probeů mezi vašimi vzorky a referencí
common_probes <- intersect(rownames(beta_m), rownames(ref_beta))
beta_all <- cbind(ref_beta[common_probes, , drop=FALSE], beta_m[common_probes, , drop=FALSE])

# 4) Převod beta -> M hodnoty (stabilizovat pomocí eps)
beta_to_M <- function(b, eps = 1e-6) log2((pmax(eps, b)) / pmax(eps, 1 - b))
M_all <- beta_to_M(beta_all)

# 5) Normalizace mezi „arrays“ (kvantily) na M scale
M_norm <- normalizeBetweenArrays(M_all, method = "quantile")

# 6) Vytvořit PCNSL centroid (pokud máte ref_labels)
pcnsl_idx <- which(ref_labels == "PCNSL")
pcnsl_centroid <- rowMeans(M_norm[, pcnsl_idx, drop=FALSE], na.rm = TRUE)

# 7) Spočítat Pearson korelaci každého vzorku (vaše 5 vzorků jsou v colnames beta_m)
sample_names <- colnames(beta_m)
similarity <- map_dbl(sample_names, ~ cor(pcnsl_centroid, M_norm[, .x], use = "pairwise.complete.obs"))
names(similarity) <- sample_names
similarity  # vektor korelací (vyšší = podobnější PCNSL)

# 8) Jednoduchý DMP test v referenci (PCNSL vs ostatní) pomocí limma — pokud máte dost referenčních vzorků
if(sum(ref_labels == "PCNSL") >= 2) {
  design <- model.matrix(~ 0 + factor(ref_labels))
  colnames(design) <- levels(factor(ref_labels))
  fit <- lmFit(M_norm[, colnames(ref_beta)], design)
  contrast <- makeContrasts(PCNSL_vs_other = PCNSL - (rowMeans(design[, setdiff(colnames(design),"PCNSL")])), levels = design)
  # Pozn: výše contrast vyžaduje úpravu podle skutečných skupin — místo toho použijte jednoduchý dvě‑skupinový design když máte jen PCNSL vs non
  # fit2 <- contrasts.fit(fit, contrasts = your_contrast_matrix)
  # fit2 <- eBayes(fit2)
  # topTable(fit2, number = 50)
}

# 9) Rychlá kontrola kandidátních genů (seznam genů u kterých chcete kontrolovat metylaci)
genes_of_interest <- c("CXCR4","S1PR1","ITGA4","ICAM1","VCAM1","MMP9","CLDN5","CD274","HLA-A","TJP1")
# Pokud máte manifest mapping probe->gene (manifest dataframe manifest with columns probe and gene), získejte jejich probe set a prohlédněte hodnoty:
# manifest %>% filter(gene %in% genes_of_interest) %>% pull(probe) -> probes_of_interest
# M_norm[probes_of_interest, sample_names, drop=FALSE]

# Výstupy: similarity (k PCNSL signatuře) a DMP výsledky z limma (pokud máte reference).