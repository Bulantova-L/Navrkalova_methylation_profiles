library(tidyverse)
setwd("methylation_profiles")

matched_data <- read_tsv("matched_cpgs.bed", col_names = FALSE)

# Assuming:
# Col 11 = Nanopore methylation percentage (0-100)
# Col 22 = Illumina Probe ID (cgXXXXXXXX)
patient_probes <- matched_data %>%
  select(probe_id = X22, nano_meth = X11) %>%
  mutate(nano_beta = nano_meth / 100) %>%
  group_by(probe_id) %>%
  summarise(patient_beta = mean(nano_beta, na.rm = TRUE))

# For multi-sample: bind rows from each sample then pivot_wider
# final_matrix <- long_df |> pivot_wider(names_from = sample, values_from = mean_meth)





# Load your reference matrix (Rows = Probes, Columns = Samples)
ref_matrix <- read.delim("reference_epic_matrix.txt", row.names = 1)

# Merge your patient data with the reference database
merged_data <- ref_matrix %>%
  rownames_to_column("probe_id") %>%
  inner_join(patient_probes, by = "probe_id") %>%
  column_to_rownames("probe_id")

# Now 'merged_data' contains your patient as a new column alongside the references!

# Load the manifest to get the Probe-to-Gene mapping
manifest <- read_tsv("EPICv2.hg38.manifest.tsv.gz", comment = "#") %>%
  select(probe_id = id, gene_name = gene_symbol) # Adjust column names based on actual file headers

target_genes <- c("CXCR4", "CXCR5", "CXCL12", "CXCL13", "CCR7", "ITGA4", "MMP2", "MMP9", "CD274", "MYD88", "CD79B")

# Get a list of probe IDs that belong to your target genes
target_probes <- manifest %>%
  filter(gene_name %in% target_genes) %>%
  pull(probe_id)

# Filter your merged matrix down to ONLY these biologically relevant probes
evaluation_matrix <- merged_data[rownames(merged_data) %in% target_probes, ]

# Separate your patient from the reference groups
patient_vector <- evaluation_matrix$patient_beta
ref_only <- evaluation_matrix[, colnames(evaluation_matrix) != "patient_beta"]

# Calculate Spearman correlation between patient and all reference samples
correlations <- cor(patient_vector, ref_only, method = "spearman", use = "complete.obs")

# Look at the distribution: Is the correlation coefficient higher against 
# the known PCNSL samples or the systemic DLBCL samples?
library(ggfortify)

# Transpose matrix so samples are rows and probes are columns
pca_data <- t(na.omit(evaluation_matrix))

# Run PCA
pca_results <- prcomp(pca_data, scale. = TRUE)

# Create a metadata data frame mapping sample names to their groups 
# (e.g., "PCNSL", "Systemic_DLBCL", "Your_Patient")
metadata <- data.frame(
  Sample = rownames(pca_data),
  Group = c(rep("Systemic_DLBCL", 20), rep("PCNSL", 20), "Your_Patient") # Example groups
)

# Plot it
autoplot(pca_results, data = metadata, colour = 'Group', size = 3) +
  theme_minimal() +
  labs(title = "Lymphoma CNS-Homing Methylation Space")