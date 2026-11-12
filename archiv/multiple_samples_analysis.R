# language: r
# Potřebné balíčky: tidyverse, limma
library(tidyverse)
library(limma)
library(dplyr)
# 1. Load all epic_mapped.bed files
file_paths <- list.files(
  path = "results",
  pattern = "_epic_mapped\\.bed$",
  recursive = TRUE,
  full.names = TRUE
)

# Read and combine all files
meth_data <- lapply(file_paths, function(f) {
  sample_name <- sub(".*/(.+?)_epic_mapped\\.bed", "\\1", f)
  read_tsv(f, col_names = c("chrom", "start", "end", "probe_id", "methylation"),
           show_col_types = FALSE) %>%
    mutate(sample = sample_name)
}) %>% bind_rows()

# 2. Pivot to wide format (probes as rows, samples as columns)
meth_wide <- meth_data %>%
  mutate(methylation = na_if(methylation, ".")) %>%
  mutate(methylation = as.numeric(methylation)) %>%
  select(probe_id, sample, methylation) %>%
  pivot_wider(names_from = sample, values_from = methylation)

# 3. Basic summary statistics
summary(meth_wide %>% select(-probe_id))

# 4. Check data dimensions and missing values
dim(meth_wide)
colSums(is.na(meth_wide))


# 2. PCA Analysis
# Remove probes with missing values and prepare matrix
meth_matrix <- meth_wide %>%
  select(-probe_id) %>%
  as.matrix()

rownames(meth_matrix) <- meth_wide$probe_id

# Remove rows with any NA
meth_complete <- meth_matrix[complete.cases(meth_matrix), ]

# Remove probes with zero or near-zero variance
probe_variance <- apply(meth_complete, 1, var, na.rm = TRUE)
meth_filtered <- meth_complete[probe_variance > 0, ]

# Check how many probes were removed
cat("Probes removed:", nrow(meth_complete) - nrow(meth_filtered), "\n")
cat("Probes remaining:", nrow(meth_filtered), "\n")

# Now run PCA on filtered data
pca_result <- prcomp(t(meth_filtered), scale. = TRUE, center = TRUE)

# PCA plot
pca_df <- data.frame(
  sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5) +
  labs(title = "PCA of Methylation Profiles",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)")) +
  theme_minimal()












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