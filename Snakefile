# language: python
import os
import pandas as pd
import subprocess

#snakemake --use-conda --cores 28

# Načtení seznamu vzorků z TSV souboru vytvořeného v R
samples_df = pd.read_csv("samples_list.tsv", sep="\t", header=None, names=["sample", "path"])
SAMPLES = samples_df["sample"].tolist()
BAM_MAP = dict(zip(samples_df["sample"], samples_df["path"]))

MANIFEST = "epic_array_probes_hg38.clean.sorted.bed"

rule all:
    input:
        expand("results/{sample}/{sample}_epic_mapped.bed", sample=SAMPLES)

rule ensure_sorted_bam:
    input:
        raw_bam = lambda wc: BAM_MAP[wc.sample]
    output:
        bam = "results/{sample}/{sample}.sorted.bam",
        bai = "results/{sample}/{sample}.sorted.bam.bai"
    threads: 8
    conda: "envs/samtools.yaml"
    shell:
        """
        # Převedení na absolutní cesty
        inp=$(readlink -f {input.raw_bam})
        out_bam=$(readlink -f {output.bam})
        out_bai=$(readlink -f {output.bai})
        
        mkdir -p $(dirname $out_bam)

        # Kontrola, zda je soubor seřazený (přes hlavičku BAMu)
        if samtools view -H "$inp" | grep -q "SO:coordinate"; then
            # UŽ JE SEŘAZENO -> Vytvoříme symlink
            rm -f "$out_bam"
            ln -s "$inp" "$out_bam"
            
            # Hledání existujícího indexu
            rm -f "$out_bai"
            if [ -f "$inp.bai" ]; then
                ln -s "$inp.bai" "$out_bai"
            elif [ -f "${{inp%.bam}}.bai" ]; then
                ln -s "${{inp%.bam}}.bai" "$out_bai"
            else
                samtools index "$out_bam"
            fi
        else
            # NENÍ SEŘAZENO -> setřídíme a zaindexujeme
            samtools sort -@ {threads} -o "$out_bam" "$inp"
            samtools index "$out_bam"
        fi
        """

rule download_manifest:
    output:
        MANIFEST
    shell:
        """
        wget -qO- https://raw.githubusercontent.com/eleozzandrini/EPIC_manifest_hg38/master/EPIC_hg38_all.bed | \
        grep "^chr" | sort -k1,1 -k2,2n > {output}
        """

rule modkit_pileup:
    input:
        bam = "results/{sample}/{sample}.sorted.bam" # Modkit teď vždy bere zajištěný BAM/symlink
    output:
        bed = "results/{sample}/{sample}_modkit.bed",
        log = "results/{sample}/{sample}_modkit.log"
    threads: 28
    conda: "envs/modkit.yaml"
    shell:
        """
        modkit pileup {input.bam} {output.bed} \
            --filter-threshold C:0.8 --filter-threshold A:0.9 \
            --threads {threads} \
            --log-filepath {output.log}
        """

rule sort_modkit_bed:
    input:
        "results/{sample}/{sample}_modkit.bed"
    output:
        "results/{sample}/{sample}_modkit_sorted.bed"
    threads: 4
    shell:
        """
        # Lexikografické třídění stávajícího modkit outputu
        LC_ALL=C sort -k1,1 -k2,2n {input} > {output}
        """

rule bedtools_map:
    input:
        manifest = MANIFEST,
        # Nyní bere jako vstup ten nově setříděný soubor:
        modkit_bed = "results/{sample}/{sample}_modkit_sorted.bed"
    output:
        mapped = "results/{sample}/{sample}_epic_mapped.bed"
    shell:
        """
        bedtools map -a {input.manifest} \
            -b {input.modkit_bed} \
            -c 11 -o mean > {output.mapped}
        """