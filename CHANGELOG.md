# Changelog
## 1.0.0 (2026-10-28)
- First release of MutaXplorer (developed from the Trivanta prototype).
- Inputs: DNA FASTQ / BAM, RNA FASTQ; hg19 / hg38; WGS / WES / panel presets.
- SNV/indel: Mutect2 (tumour-only, tumour-normal), HaplotypeCaller; UMI consensus (fgbio); BQSR.
- New: CNV (CNVkit, gene-level calls), SV (Manta somatic + germline), DNA fusions from SV breakpoints,
  RNA fusions (STAR + Arriba) with DNA/RNA evidence integration, SBS96 signatures with NNLS fitting,
  LOH-robust sample-identity check, oncoplot with CNV/fusion markers.
- TMB, MSI (MSIsensor2), VEP annotation, MultiQC, cohort tables, Excel and HTML report.
- Simulated test suite (4 scenarios) and recall benchmark; CI workflow.
