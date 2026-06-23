# Agentic workflow — methylation_profiles

Purpose
- Short guide for automated/agentic interactions with the methylation_profiles folder: quick run commands, key files, and minimal operational notes.

Contents (key files)
- cpgIslandExt.txt — CpG island annotations
- EPICv2.hg38.manifest.tsv.gz — probe manifest (gzipped)
- epic_array_probes_hg38.clean.bed, .clean.sorted.bed — cleaned probe BED (sorted)
- epic_array_probes_hg38.bed — original probe BED
- hg38.chrom.sizes, hg38_cpg_islands.bed — reference resources
- matched_cpgs.bed — matched CpG positions used by analyses
- methylation/ — directory with methylation input and results
- methylation_profile_test.sh — smoke/test script
- extract_all_promoters.R — promoter extraction helper
- statistical_analysis_cpgs.R — core analysis script

Quick start
- Change to repository folder and run the smoke tests:

```bash
# Shell
cd methylation_profiles
bash methylation_profile_test.sh
```

- Run the main R analysis (non-interactive):

```bash
# R
Rscript --vanilla statistical_analysis_cpgs.R
```

- Regenerate sorted BED if needed:

```bash
# Shell
sort -k1,1 -k2,2n epic_array_probes_hg38.clean.bed > epic_array_probes_hg38.clean.sorted.bed
```

Agentic steps for automation
1. Ensure environment: R >= 4.0 and required R packages (see scripts for library calls). Install bedtools/tabix if BED/VCF utilities are used.
2. Run lint/smoke: methylation_profile_test.sh. Fail early on missing inputs.
3. Run extract_all_promoters.R when promoter regions need refreshing.
4. Run statistical_analysis_cpgs.R to produce results; have outputs directed to methylation/.
5. Archive outputs and manifest versions (keep EPICv2.hg38.manifest.tsv.gz alongside results).

Notes
- Scripts assume working directory = methylation_profiles. Use absolute paths in automation where practical.
- Check script headers for package requirements and any hard-coded paths before running in a CI/agent environment.

Actions I can take
- Add a CI job, create a detailed README, or generate a Dockerfile for reproducible runs on CI or agents.
