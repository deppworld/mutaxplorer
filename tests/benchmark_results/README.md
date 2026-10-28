# Simulation benchmark (manuscript Figure 2 and benchmark table)

Reproduce (~1 h on 2 CPUs):
```bash
python3 tests/benchmark_simulation.py make bench --depth 500
nextflow run . --input bench/samplesheet_paired.csv    --fasta bench/genome.fa --targets bench/targets.bed \
     --assay panel --analysis somatic --germline_resource bench/gnomad.vcf.gz --outdir bench/results_paired
nextflow run . --input bench/samplesheet_tumoronly.csv --fasta bench/genome.fa --targets bench/targets.bed \
     --assay panel --analysis somatic --germline_resource bench/gnomad.vcf.gz --outdir bench/results_tumoronly -resume
python3 tests/benchmark_simulation.py eval bench bench/results_paired bench/results_tumoronly bench/benchmark
```
`benchmark.tsv` = per-run recall / germline leakage / precision; `benchmark.png` = summary figure.
Run with MutaXplorer 1.0.0, GATK 4.6.2.0, bwa-mem2, fastp, samtools/bcftools 1.24 (see pipeline_info).
