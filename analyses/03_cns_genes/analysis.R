library(biomaRt)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# CNS geny
cns_genes <- c(
  "CXCR4", "CXCL12", "S1PR1", "S1PR2",  # Migrace
  "MMP9", "CLDN5", "TJP1",               # BBB
  "CD274", "HLA-A", "B2M",               # Imunitní únik
  "MYD88", "CD79B", "PIM1", "PRDM1"      # Onkogeny
)

options(timeout = 300)                # 300 s
# Připojit k Ensembl (hg38 - odpovídá vašemu Snakefile)
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl",
                      version = 109)  # GRCh38

# Získat genomické koordináty
gene_coords <- getBM(
  attributes = c('hgnc_symbol', 'chromosome_name', 
                 'transcription_start_site', 'strand',
                 'start_position', 'end_position'),
  filters = 'hgnc_symbol',
  values = cns_genes,
  mart = ensembl
)

# Extrahovat CpG v promoterech
extract_cns_cpgs <- function(all_data, gene_regions) {
  cns_cpgs <- all_data |>
    separate(probe_id, into = c("chr", "pos"), sep = ":", remove = FALSE) |>
    mutate(pos = as.numeric(pos))
  
  # Filter pro každý gen
  results <- list()
  
  for(i in 1:nrow(gene_regions)) {
    gene <- gene_regions[i, ]
    
    gene_cpgs <- cns_cpgs |>
      filter(
        chr == gene$chr,
        pos >= gene$promoter_start,
        pos <= gene$promoter_end
      )
    
    if(nrow(gene_cpgs) > 0) {
      gene_cpgs$gene <- gene$hgnc_symbol
      results[[gene$hgnc_symbol]] <- gene_cpgs
    }
  }
  
  bind_rows(results)
}

cns_methylation <- extract_cns_cpgs(all_samples, gene_coords)


# Vytvořit promoter regiony (-2kb/+500bp od TSS)
cns_regions <- gene_coords |>
  filter(chromosome_name %in% c(1:22, "X", "Y")) |>
  mutate(
    chr = paste0("chr", chromosome_name),
    # Promoter
    promoter_start = ifelse(strand == 1, 
                           transcription_start_site - 2000,
                           transcription_start_site - 500),
    promoter_end = ifelse(strand == 1,
                         transcription_start_site + 500,
                         transcription_start_site + 2000),
    # Celý gen (pro analýzu genového těla)
    gene_start = start_position,
    gene_end = end_position
  )

# Připravit BED soubor - PROMOTERY
cns_promoters_bed <- cns_regions |>
  dplyr::select(chr, promoter_start, promoter_end, hgnc_symbol) |>
  arrange(chr, promoter_start)

# Zapsat BED soubor
fwrite(cns_promoters_bed, 
       "data/cns_promoters_hg38.bed",   # vstup pro Snakefile (CNS_MANIFEST)
       sep = "\t", 
       col.names = FALSE)

cat("Vytvořeno:", nrow(cns_promoters_bed), "CNS promoter regionů\n")

# Také celé genové regiony (pro komplexnější analýzu)
cns_genes_bed <- cns_regions |>
  dplyr::select(chr, gene_start, gene_end, hgnc_symbol) |>
  arrange(chr, gene_start)

fwrite(cns_genes_bed, 
       "data/cns_genes_hg38.bed",
       sep = "\t", 
       col.names = FALSE)

cat("Vytvořeno:", nrow(cns_genes_bed), "CNS geny\n")



# Fáze 1: Stažení referenčních dat z GEO
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("minfi")

library(GEOquery)
library(minfi) # pro práci s EPIC array daty
library(tidyverse)
library(data.table)

# PCNSL dataset GSE25297, GSE28952
# DLBCL GSE145043

# Příklad: PCNSL dataset
gse <- getGEO("GSE25297")
ref_pheno <- pData(phenoData(gse[[1]]))
# Stáhnout raw IDAT soubory nebo již zpracovaná beta values


# gse <- getGEO("GSE145043", GSEMatrix = TRUE)
# metadata <- pData(phenoData(gse[[1]]))

# # Filter only the samples that are eRRBS (Methylation)
# # Based on the GEO record, they are named "DLBCL tissue eRRBS"
# methylation_samples <- metadata %>%
#   filter(str_detect(`characteristics_ch1`, "eRRBS") | str_detect(title, "eRRBS")) %>%
#   select(geo_accession, title, characteristics_ch1)

# head(methylation_samples)

# # Download supplementary file for a specific GSM accession
# gsm_files <- getGEOSuppFiles("GSM4304988")

# # Find the path of the downloaded file
# file_path <- rownames(gsm_files)[1]

# # Read the file using data.table (fastest for large genomics text files)
# # Bismark/eRRBS outputs usually contain: chr, start, end, methylation percentage, coverage
# meth_data <- fread(file_path)
# head(meth_data)

download_errbs_samples <- function(gsm_ids, output_dir = "data/geo_data") {
  dir.create(output_dir, showWarnings = FALSE)
  
  sample_list <- list()
  
  for (gsm_id in gsm_ids) {
    message("Downloading ", gsm_id)
    
    tryCatch({
      gsm_files <- getGEOSuppFiles(gsm_id, baseDir = output_dir)
      file_path <- rownames(gsm_files)[1]
      
      # Načíst data
      meth <- fread(file_path)
      
      # Standardizovat formát: chr:position jako ID
      meth[, probe_id := paste0(chr, ":", base)]
      meth[, beta := freqC / 100]  # Konverze na 0-1 škálu
      
      # Filtrovat pokrytí >= 10
      meth <- meth[coverage >= 10]
      
      sample_list[[gsm_id]] <- meth[, .(probe_id, beta)]
      
    }, error = function(e) {
      message("Failed: ", gsm_id, " - ", e$message)
    })
    
    Sys.sleep(2)  # Rate limiting
  }
  
  sample_list
}

test_samples <- head(methylation_samples$geo_accession, 10)
dlbcl_data <- download_errbs_samples(test_samples)

# Spojit všechny vzorky do wide formátu
library(data.table)

# Získat společné CpG pozice napříč všemi vzorky
merge_samples <- function(sample_list) {
  # Start s prvním vzorkem
  merged <- sample_list[[1]]
  setnames(merged, "beta", names(sample_list)[1])
  
  # Postupně přidat další vzorky
  for(i in 2:length(sample_list)) {
    sample_name <- names(sample_list)[i]
    temp <- sample_list[[i]]
    setnames(temp, "beta", sample_name)
    merged <- merge(merged, temp, by = "probe_id", all = FALSE)  # inner join
  }
  
  merged
}

dlbcl_matrix <- merge_samples(dlbcl_data)

# Kolik společných CpG pozic?
cat("Počet společných CpG:", nrow(dlbcl_matrix), "\n")
cat("Počet vzorků:", ncol(dlbcl_matrix) - 1, "\n")


# 1. Načíst plná nanopore data (ne epic_mapped.bed, ale modkit.bed)
library(data.table)

# Cesta k vašim původním modkit výstupům
nanopore_files <- list.files(
  path = "results", 
  pattern = "_modkit_sorted.bed$", 
  recursive = TRUE, 
  full.names = TRUE
)

# Načíst vaše nanopore vzorky
read_nanopore <- function(file_path) {
  sample_name <- str_extract(basename(file_path), "^[^_]+")
  
  dt <- fread(file_path, select = c(1, 2, 11))  # chr, position, % methylation
  setnames(dt, c("chr", "pos", sample_name))
  
  # Vytvořit probe_id ve formátu chr:pos
  dt[, probe_id := paste0(chr, ":", pos)]
  dt[, c("probe_id", sample_name), with = FALSE]
}

# Načíst všechny vaše vzorky
your_samples <- lapply(nanopore_files, read_nanopore)

# Spojit do jedné matice
your_matrix <- Reduce(function(x, y) merge(x, y, by = "probe_id", all = FALSE), 
                      your_samples)

cat("Vaše nanopore data - CpG:", nrow(your_matrix), "\n")



# Vaše data už máte v meth_wide - spojit s DLBCL daty
# Najít společné probe_id

# Převést vaše chr:pos na stejný formát
# (předpokládám že meth_wide má chr1:12345 formát)

# Spojit datasety
all_samples <- merge(
  dlbcl_matrix, 
  your_data, 
  by = "probe_id", 
  all = FALSE  # jen společné CpG
)

cat("Finální počet CpG:", nrow(all_samples), "\n")
cat("Celkem vzorků:", ncol(all_samples) - 1, "\n")