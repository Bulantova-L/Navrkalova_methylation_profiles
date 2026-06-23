# Run this inside R
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "rtracklayer"))

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)

# 1. Get all gene models for hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 2. Define promoter region (TSS -2000 to +200)
promoters_gr <- promoters(txdb, upstream = 2000, downstream = 200)

# 3. Swap the cryptic Entrez IDs for standard Gene Symbols (e.g., "CXCR4")
gene_ids <- names(promoters_gr)
symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID")
promoters_gr$gene_symbol <- symbols

# 4. Drop genes that don't have a valid symbol mapping or are on weird alternative scaffolds
promoters_gr <- promoters_gr[!is.na(promoters_gr$gene_symbol)]
promoters_gr <- keepStandardChromosomes(promoters_gr, pruning.mode = "coarse")

# 5. Convert to data frame to save cleanly as a proper BED file
promoters_df <- as.data.frame(promoters_gr)
bed_output <- promoters_df[, c("seqnames", "start", "end", "gene_symbol")]

# Adjust R's 1-based coordinates to BED format's 0-based coordinates
bed_output$start <- bed_output$start - 1

# Write out the file
write.table(bed_output, "hg38_promoters.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)