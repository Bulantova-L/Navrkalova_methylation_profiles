import pandas as pd
import glob
import os

# --- Configuration ---
MIN_COVERAGE = 5
WINDOW_SIZE = 100000
# Replace with your actual modkit bedMethyl files
BEDMETHYL_FILES = glob.glob("results/*/methylation/*_modkit_sorted.bed") 

def process_sample(file_path):
    """Parses a modkit bedMethyl file and returns mean methylation per window."""
    
    # Modkit bedMethyl relevant columns:
    # 0: chrom, 1: start, 9: valid_coverage, 10: fraction_modified (0-100%)
    df = pd.read_csv(
        file_path, sep='\t', header=None, 
        usecols=[0, 1, 9, 10], 
        names=['chrom', 'start', 'coverage', 'fraction_modified']
    )
    
    # 1. Filter out low confidence sites
    df = df[df['coverage'] >= MIN_COVERAGE].copy()
    
    # 2. Convert percentage (0-100) to a fraction (0-1)
    df['ratio'] = df['fraction_modified'] / 100.0
    
    # 3. Assign each CpG site to a 100kb window block
    df['window_start'] = (df['start'] // WINDOW_SIZE) * WINDOW_SIZE
    df['window_id'] = df['chrom'] + '_' + df['window_start'].astype(str)
    
    # 4. Calculate the depth-agnostic mean ratio per window
    window_means = df.groupby('window_id')['ratio'].mean()
    
    # Label the data with the sample name
    sample_name = os.path.basename(file_path).replace('.bed', '')
    window_means.name = sample_name
    
    return window_means

# --- Main Execution ---
print("Processing modkit samples...")
sample_series = [process_sample(f) for f in BEDMETHYL_FILES]

# Combine all sample series into one large matrix (columns = samples, rows = windows)
pca_matrix = pd.concat(sample_series, axis=1)

# 5. Drop windows that lack sufficient coverage in one or more samples
pca_matrix_clean = pca_matrix.dropna()

# Save the final matrix
pca_matrix_clean.to_csv("analyses/01_cpg_counts_windows/output/pca_methylation_matrix.tsv", sep='\t', index_label='window_id')
print(f"Matrix saved! Final dimensions for PCA: {pca_matrix_clean.shape}")