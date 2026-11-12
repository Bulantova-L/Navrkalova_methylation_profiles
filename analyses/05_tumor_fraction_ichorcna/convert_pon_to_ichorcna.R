# Prevod externiho PoN (data.frame chr/start/end/Median) do formatu, ktery cte ichorCNA
#
# ichorCNA (runIchorCNA.R -> normalizeByPanelOrMatchedNormal) ocekava:
#   - GRanges (ne data.frame)
#   - biny presne shodne s biny vzorku (findOverlaps type = "equal"): 1-based, start = k*1e6 + 1, end = (k+1)*1e6
#   - Median v log2 skale (log2 korigovane hloubky, ~0 pro normalni bin)
# ONT_PoN_1control_1Mb.rds ma end = start + 1e6 (o 1 bp delsi) a Median v linearni skale (~1)
#   -> zadny bin se neshodne, ichorCNA dostane 0 binu a pada v runEM ("rho[, 1]: subscript out of bounds")
#
# Spusteni (pracovni adresar = methylation_profiles, conda activate ichorcna):
#   Rscript analyses/05_tumor_fraction_ichorcna/convert_pon_to_ichorcna.R \
#     analyses/05_tumor_fraction_ichorcna/pon/ONT_PoN_1control_1Mb.rds analyses/05_tumor_fraction_ichorcna/pon/ONT_PoN_1control_1Mb.GR.rds

suppressPackageStartupMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Pouziti: Rscript convert_pon_to_ichorcna.R <vstup.rds> <vystup.rds>")
MIN_RATIO <- 0.1   # biny s linearnim pomerem < 0.1 v kontrole = artefakt (nulove pokryti, chrY u zeny) -> NA

x <- as.data.frame(readRDS(args[1]))
stopifnot(all(c("chr", "start", "end", "Median") %in% names(x)))

# 1) Biny: end o 1 bp kratsi, pokud jsou o 1 bp delsi nez krok
w <- x$end - x$start
if (all(w %% 1000 == 0)) x$end <- x$end - 1L
message("Sirka binu po uprave: ", paste(unique(x$end - x$start + 1), collapse = ", "))

# 2) Skala: linearni pomer -> log2
med <- median(x$Median, na.rm = TRUE)
if (med > 0.5) {
  message("Median panelu = ", round(med, 3), " -> linearni skala, prevadim na log2")
  bad <- !is.na(x$Median) & x$Median < MIN_RATIO
  message("Binu s pomerem < ", MIN_RATIO, " (-> NA): ", sum(bad))
  x$Median[bad] <- NA
  x$Median <- log2(x$Median)
} else {
  message("Median panelu = ", round(med, 3), " -> uz log2, beze zmeny")
}

gr <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
saveRDS(gr, args[2])
print(gr)
print(summary(gr$Median))
