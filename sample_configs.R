# language: r
setwd("z:/share/710000-CEITEC/713000-cmm/713016-bioit/base/workspace/lucka/Navrkalova_cfDNA/methylation_profiles")
common_path <- "/mnt/data/ceitec_cfg2/710000-CEITEC/713000-cmm/713003-pospisilova/base/MedGen/ONTdata/Lympho_ONT"

# Definice cest
bam_dict <- tibble(
  sample_name = c("cfDNA_Pilot1", "cfDNA_BAR", "cfDNA_LOS", "cfDNA_MUD", "cfDNA_VN"),
  bam_path = c(
    file.path(common_path, "cfDNA_Pilot1/cfDNA_DLBCL_sup_mC_aligned_output.sorted.bam"),
    file.path(common_path, "cfDNA_DLBCL_BAR,LOS/cfDNA_BAR/bar.sorted.bam"),
    file.path(common_path, "cfDNA_DLBCL_BAR,LOS/cfDNA_LOS/los.bam"),
    file.path(common_path, "cfDNA_MUD/cfDNA_MUD_hg38.bam"),
    file.path(common_path, "cfDNA_VN/cfDNA_VN_hg38.bam")
  )
)

# Uložení pro Bash
write_tsv(bam_dict, "samples_list.tsv", col_names = FALSE)
