# Test for methylation profile generation using a BAM file as input
conda activate modkit

# 1. Kontrola, zda byl zadán argument (jméno vzorku)
TARGET_SAMPLE=$1
LIST_FILE="samples_list.tsv"
MANIFEST="epic_array_probes_hg38.clean.sorted.bed"
OUTPUT_BASE="results"


if [ -z "$TARGET_SAMPLE" ]; then
    echo "Nebyl zadán konkrétní vzorek, zpracovávám vše ze $LIST_FILE"
    # Vybere všechny řádky
    DATA_TO_PROCESS=$(cat "$LIST_FILE")
else
    echo "Hledám vzorek: $TARGET_SAMPLE"
    # Vyhledá v souboru pouze řádek pro daný vzorek (přesná shoda na 1. sloupci)
    DATA_TO_PROCESS=$(awk -v target="$TARGET_SAMPLE" '$1 == target {print $0}' "$LIST_FILE")
    
    if [ -z "$DATA_TO_PROCESS" ]; then
        echo "Chyba: Vzorek '$TARGET_SAMPLE' nebyl v $LIST_FILE nalezen."
        exit 1
    fi
fi

# --- KONTROLA MANIFESTU ---
if [ ! -f "$MANIFEST" ]; then
    echo "Manifest '$MANIFEST' nebyl nalezen v aktuálním adresáři."
    ## Step 1: Prepare the EPIC array manifest for hg38
    # wget https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationV1/main/Anno/EPICv2/EPICv2.hg38.manifest.tsv.gz
    gzip -dc EPICv2.hg38.manifest.tsv.gz | awk 'BEGIN{OFS="\t"} NR>1 && $1!~/^#/ {print $1,$2,$3,$9}' > epic_array_probes_hg38.bed
    awk 'BEGIN{OFS="\t"} $1!~/^NA$/ && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $2>=0 && $3>=0 {print $1,$2,$3,$4}' epic_array_probes_hg38.bed > epic_array_probes_hg38.clean.bed
    sort -k1,1 -k2,2n epic_array_probes_hg38.clean.bed > epic_array_probes_hg38.clean.sorted.bed
else
    echo "Používám existující manifest: $MANIFEST"
fi

# --- LOOP PŘES VZORKY ---
echo "$DATA_TO_PROCESS" | while IFS=$'\t' read -r name path; do
    echo "--------------------------------------------------------"
    echo "VZOREK: $name"
    echo "BAM: $path"
    
    SAMPLE_DIR="$OUTPUT_BASE/$name"
    mkdir -p "$SAMPLE_DIR"

    output_bed="$SAMPLE_DIR/methylation/${name}_modkit.bed"
    # 1. Modkit pileup (volitelně: přidejte logování do souboru)
    echo "[1/2] Spouštím modkit pileup..."
    modkit pileup "$path" $output_bed \
        --filter-threshold C:0.8 --filter-threshold A:0.9 \
        --threads 25 \
        --log-filepath "$SAMPLE_DIR/methylation/${name}_modkit.log"

    # 2. Agregace na EPIC sondy přes bedtools
    echo "[2/2] Mapování na EPIC manifest..."
    bedtools map -a "$MANIFEST" \
        -b $output_bed \
        -c 11 -o mean > "$SAMPLE_DIR/${name}_epic_mapped.bed"

    echo "HOTOVO: Výsledky uloženy v $SAMPLE_DIR"
done




# # parallel jobs
# N=8
# tmp=$(mktemp -d /tmp/bedparts.XXXX)

# # Method A — split by chromosome (recommended)
# awk 'BEGIN{OFS="\t"}{print > "'"$tmp"'/A."$1".bed"}' methylation/${sample_name}/${sample_name}_modkit.bed
# ls "$tmp"/A.*.bed | parallel -j "$N" 'bedtools intersect -sorted -a {} -b epic_array_probes_hg38.clean.sorted.bed -wa -wb > {}.int'

# Method B — (alternative) split by line-count for even work distribution
# split -l 1000000 output.sorted.bed "$tmp"/chunk.
# ls "$tmp"/chunk.* | parallel -j "$N" 'bedtools intersect -sorted -a {} -b epic_array_probes_hg38.clean.sorted.bed -wa -wb > {}.int'

# # Merge results and cleanup
# cat "$tmp"/*.int > matched_cpgs.bed
# rm -r "$tmp"


# # compute mean methylation per EPIC probe (requires sorted files)
# bedtools map -a epic_array_probes_hg38.clean.sorted.bed \
#               -b methylation/${sample_name}/${sample_name}_modkit.bed \
#               -c 11 -o mean > epic_probe_meth_mean.bed

# # Downloading CpG Islands from UCSC
# wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz
# gunzip cpgIslandExt.txt.gz

# # Convert UCSC format to standard BED format (chrom, start, end, name)
# awk 'BEGIN {OFS="\t"} {print $2, $3, $4, $5}' cpgIslandExt.txt > hg38_cpg_islands.bed

# ## Aggregate the Methylation Data
# # -a: Your functional windows (Promoters)
# # -b: Your single-CpG nanopore data
# # -c 11: Apply the operation on column 11 (methylation frequency)
# # -o mean: Calculate the average of all overlapping points

# # Now run the R script for extracting promoter regions
# Rscript extract_all_promoters.R

# bedtools map -a hg38_promoters.bed -b ${output_bed} -c 11 -o mean > nanopore_aggregated_promoters.bed