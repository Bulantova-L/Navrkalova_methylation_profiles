# language: python
import os
import pandas as pd
import subprocess

METHYLATION = {
    "5mC": "m",
    "6mA": "a",
    "5hmC": "h"}

basecaller_location = "/tmp/dorado-0.5.3-linux-x64/bin/dorado"

#snakemake --use-conda --cores 28
# snakemake --use-conda --cores 28 2>&1 | tee snakemake.log


# Načtení seznamu vzorků z TSV souboru vytvořeného v R
samples_df = pd.read_csv("data/samples_list.tsv", sep="\t", header=None, names=["sample", "path"])
SAMPLES = samples_df["sample"].tolist()
BAM_MAP = dict(zip(samples_df["sample"], samples_df["path"]))

MANIFEST = "data/epic_array_probes_hg38.clean.sorted.bed"
CNS_MANIFEST = "data/cns_promoters_hg38.bed"
# hg38 CpG loci from GSE54200 (Mann-Whitney probes), lifted over + sorted by
# analyses/07_gse54200_dlbcl/DLBCL_metylace.R + DLBCL_metylace.sh (formerly ../DLBCL_metylace)
# If the pipeline is ever re-run, use `snakemake --touch` first so the moved file doesn't trigger a re-run.
GSE54200_MANIFEST = "analyses/07_gse54200_dlbcl/reference/output_hg38_sorted.bed"

rule all:
    input:
        #expand("results/{sample}/methylation/{sample}_epic_mapped.bed", sample=["cfDNA_VN", "cfDNA_MUD"]),
        #expand("results/{sample}/methylation/{sample}_cns_mapped.bed", sample=SAMPLES),
        #expand("results/{sample}/methylation/{sample}_gse54200_mapped.bed", sample=SAMPLES),
        expand("results/{sample}/methylation/{sample}_filtered_{mod_name}.bed", sample=SAMPLES, mod_name=METHYLATION.keys()),
        expand('results/{sample}/summary/{sample}_summary.tsv', sample=SAMPLES),
        #expand('results/{sample}/qc_reports/{sample}_pycoQC.html', sample="cfDNA_VN")
        expand('results/{sample}/qc_reports/nanoplot_from_bam/{sample}_NanoPlot-report.html', sample=SAMPLES),
        expand('results/{sample}/methylation/{sample}.mosdepth.summary.txt', sample=SAMPLES),
        'results/qc_reports/all_samples/mosdepth_coverage_dist.html',
        expand("results/{sample}/methylation/{sample}_gse54200_mapped.bed", sample=SAMPLES)



rule ensure_sorted_bam:
    input:
        raw_bam = lambda wc: BAM_MAP[wc.sample]
    output:
        bam = "results/{sample}/methylation/{sample}.sorted.bam",
        bai = "results/{sample}/methylation/{sample}.sorted.bam.bai"
    threads: 8
    resources:
        mem_gb = 16
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
            # OPTIMIZATION: -m 2G restricts memory per thread. 8 threads * 2G = 16GB limit.
            samtools sort -@ {threads} -m 2G -o "$out_bam" "$inp"
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
        bam = "results/{sample}/methylation/{sample}.sorted.bam" # Modkit teď vždy bere zajištěný BAM/symlink
    output:
        bed = "results/{sample}/methylation/{sample}_modkit.bed",
        log = "results/{sample}/methylation/{sample}_modkit.log"
    threads: 12
    resources:
        mem_gb = 40
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
        "results/{sample}/methylation/{sample}_modkit.bed"
    output:
        "results/{sample}/methylation/{sample}_modkit_sorted.bed"
    threads: 12
    resources:
        # Allocate total RAM for the rule (e.g., 32 GB)
        mem_mb=32000
    shadow: "shallow"  # Redirects unique execution to a disk-backed temp folder
    shell:
        """
        # Calculate a safe buffer for GNU sort (approx 85% of total rule memory)
        # to leave breathing room for the OS and prevent OOM (Exit code 137)
        SORT_MEM_M=$(( {resources.mem_mb} * 85 / 100 ))

        # Create a local tmp directory inside the shadow work dir
        mkdir -p ./sort_tmp

        # Run memory-constrained, parallelized, disk-backed sort
        LC_ALL=C sort \
            -S "${{SORT_MEM_M}}M" \
            --parallel={threads} \
            -T ./sort_tmp \
            -k1,1 -k2,2n \
            {input} > {output}
        """

rule bedtools_map:
    input:
        manifest = MANIFEST,
        # Nyní bere jako vstup ten nově setříděný soubor:
        modkit_bed = "results/{sample}/methylation/{sample}_modkit_sorted.bed"
    output:
        mapped = "results/{sample}/methylation/{sample}_epic_mapped.bed"
    conda: "envs/modkit.yaml"
    shell:
        """
        bedtools map -a {input.manifest} \
            -b {input.modkit_bed} \
            -c 11 -o mean > {output.mapped}
        """

rule download_dorado_tar:
    output: '/tmp/dorado-0.5.3-linux-x64.tar.gz'
    shell:
        f"""
        wget -P /tmp https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.3-linux-x64.tar.gz
        """
rule extract_dorado:
    input: '/tmp/dorado-0.5.3-linux-x64.tar.gz'
    output: basecaller_location
    shell:
        f"""
        cd /tmp ;
        tar -xf {input}
        """

rule sequencing_summary:
    input: 
        basecaller_location = basecaller_location,
        bam = 'results/{sample}/methylation/{sample}.sorted.bam'
    output: 'results/{sample}/summary/{sample}_summary.tsv'
    shell:
        """
        {input.basecaller_location} summary {input.bam} > {output}
        """

rule sequencing_QC:
    input: 
        'results/{sample}/summary/{sample}_summary.tsv'
    output: 
        html = 'results/{sample}/qc_reports/{sample}_pycoQC.html',
        json = 'results/{sample}/qc_reports/{sample}_pycoQC.json'
    conda:
        "envs/QC.yaml" 
    shell:
        """
        pycoQC -f {input} -o {output.html} --json_out {output.json}
        """

rule nanoplot_QC:
    input:
        bam = "results/{sample}/methylation/{sample}.sorted.bam",
        bai = "results/{sample}/methylation/{sample}.sorted.bam.bai"
    output:
        html = "results/{sample}/qc_reports/nanoplot_from_bam/{sample}_NanoPlot-report.html",
        stats = "results/{sample}/qc_reports/nanoplot_from_bam/{sample}_NanoStats.txt"
    params:
        outdir = "results/{sample}/qc_reports/nanoplot_from_bam",
        prefix = "{sample}_"
    threads: 25
    conda: "envs/nanoplot.yaml"
    shell:
        """
        NanoPlot --bam {input.bam} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --threads {threads} \
            --N50 \
            --loglength \
            --huge
        """

rule filter_modifications:
    # Filter modifications from modkit output and at least 5 reads
    input:      
        bed = "results/{sample}/methylation/{sample}_modkit_sorted.bed"
    output: 
        filtered = "results/{sample}/methylation/{sample}_filtered_{mod_name}.bed"
    params: 
        mod_char = lambda wildcards: METHYLATION[wildcards.mod_name]
    shell:
        """
        awk '$4 == "{params.mod_char}" && $11 > 5 && $12 >= 5' {input.bed} >  {output.filtered}
        """

rule bedtools_cns_aggregate:
    input:
        cns_bed = CNS_MANIFEST,
        modkit_bed = "results/{sample}/methylation/{sample}_modkit_sorted.bed"
    output:
        mapped = "results/{sample}/methylation/{sample}_cns_mapped.bed"
    conda: "envs/modkit.yaml"
    shell:
        """
        # Průměrná metylace pro každý CNS gen
        bedtools map \
            -a {input.cns_bed} \
            -b {input.modkit_bed} \
            -c 11 -o mean,count > {output.mapped}
        """

rule bedtools_gse54200_map:
    input:
        manifest = GSE54200_MANIFEST,
        modkit_bed = "results/{sample}/methylation/{sample}_modkit_sorted.bed"
    output:
        mapped = "results/{sample}/methylation/{sample}_gse54200_mapped.bed"
    conda: "envs/modkit.yaml"
    shell:
        """
        # m and h rows for the same CpG share the same Nvalid_cov (col10) -- it's the
        # total coverage at that position, not per-mod-code -- so first collapse each
        # position to a single row: Nmod = Nmod_m + Nmod_h, Nvalid_cov taken once,
        # then map onto GSE54200 loci and compute the coverage-weighted %5mC+5hmC.
        bedtools map -a {input.manifest} \
            -b <(awk -F'\t' -v OFS='\t' '
                  $4 == "m" || $4 == "h" {{
                      key = $1 OFS $2
                      if (key != prevkey) {{
                          if (prevkey != "") print prevchrom, prevstart, prevend, nmod, ncov
                          nmod = 0; ncov = $10
                      }}
                      nmod += $12
                      prevkey = key; prevchrom = $1; prevstart = $2; prevend = $3
                  }}
                  END {{ if (prevkey != "") print prevchrom, prevstart, prevend, nmod, ncov }}
                ' {input.modkit_bed}) \
            -c 4,5 -o sum,sum -null 0 \
        | awk -F'\t' 'BEGIN{{OFS="\t"}} {{pct=($NF>0)?100*$(NF-1)/$NF:"NA"; print $0, pct}}' \
        > {output.mapped}
        """


rule mosdepth:
    input:
        bam = "results/{sample}/methylation/{sample}.sorted.bam",
        bai = "results/{sample}/methylation/{sample}.sorted.bam.bai"
    output:
        summary = "results/{sample}/methylation/{sample}.mosdepth.summary.txt",
        dist = "results/{sample}/methylation/{sample}.mosdepth.global.dist.txt"
    params:
        prefix = "results/{sample}/methylation/{sample}"
    threads: 10
    conda: "envs/mosdepth.yaml"
    shell:
        """
        mosdepth -t {threads} {params.prefix} {input.bam}
        """

rule download_plot_dist:
    output:
        "scripts/plot-dist.py"
    shell:
        """
        mkdir -p scripts
        wget -qO {output} https://raw.githubusercontent.com/brentp/mosdepth/master/scripts/plot-dist.py
        """

rule mosdepth_plot_dist:
    input:
        script = "scripts/plot-dist.py",
        dist = expand("results/{sample}/methylation/{sample}.mosdepth.global.dist.txt", sample=SAMPLES)
    output:
        html = "results/qc_reports/all_samples/mosdepth_coverage_dist.html"
    shell:
        """
        python3 {input.script} -o {output.html} {input.dist}
        """

rule alignment_multiqc:
    input: expand('results/{sample}/qc_reports/{sample}_NanoStats.txt.json', sample=SAMPLES)
    output: html= "results/qc_reports/all_samples/multiqc_report.html"
    params: outdir = "results/qc_reports/all_samples"
    conda:
        "envs/alignment_multiqc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        multiqc {input} -o {params.outdir}
        """
