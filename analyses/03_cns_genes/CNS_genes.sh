#!/bin/bash
# Aktivovat conda prostředí
# Spoustet z methylation_profiles/: bash analyses/03_cns_genes/CNS_genes.sh  (vystup data/cns_promoters_hg38.bed)

CNS_GENES="CXCR4,CXCL12,S1PR1,S1PR2,MMP9,CLDN5,TJP1,CD274,HLA-A,B2M,MYD88,CD79B,PIM1,PRDM1"

# Stáhnout RefSeq geny pro CNS geny
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -D hg38 -e \
    "SELECT chrom, txStart, txEnd, name2, strand 
     FROM refGene 
     WHERE name2 IN ('${CNS_GENES//,/\',\'}')
     ORDER BY chrom, txStart" | \
    tail -n +2 | \
    awk 'BEGIN{OFS="\t"}{
        # Promoter: -2kb/+500bp od TSS
        if ($5 == "+") {
            tss = $2
            start = (tss - 2000 > 0) ? tss - 2000 : 0
            end = tss + 500
        } else {
            tss = $3
            start = (tss - 500 > 0) ? tss - 500 : 0
            end = tss + 2000
        }
        print $1, start, end, $4
    }' | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - -c 4 -o distinct > data/cns_promoters_hg38.bed

echo "✓ Vytvořeno CNS promoter regionů: $(wc -l < data/cns_promoters_hg38.bed)"
head data/cns_promoters_hg38.bed