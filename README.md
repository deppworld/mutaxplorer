# MutaXplorer

**Comprehensive DNA + RNA variant exploration for WGS, WES and targeted panels.**
One open-source Nextflow pipeline from **FASTQ or BAM** on **hg19 or hg38** to an interactive cohort report:
SNVs/indels (somatic and germline), copy-number alterations, structural variants, DNA and RNA fusions,
TMB, MSI, mutational signatures, sample-identity checks and QC.

```
 DNA FASTQ ─ fastp ─ bwa-mem2 ─ markdup | fgbio UMI consensus ─ [BQSR] ─┐
 DNA BAM (DRAGEN / ICA / any aligner) ─ contig + read-group checks ───────┤──► analysis-ready BAMs
 RNA FASTQ ─ STAR ─ Arriba ─────────────────────────────── RNA fusions    │
                                                                           │
   SNV/indel   Mutect2 (tumour-only | tumour-normal), HaplotypeCaller ◄────┤
   CNV         CNVkit (pooled-normal or flat reference)  ◄─────────────────┤
   SV          Manta (somatic + germline) ─► DNA fusions (gene-annotated) ◄┤
   MSI         MSIsensor2        QC: mosdepth, samtools, MultiQC ◄─────────┤
   Identity    common-SNP genotypes, LOH-robust pairwise check ◄───────────┘
        │
   bcftools norm ─► Ensembl VEP ─► cohort tables · TMB · MSI status · gene-level CNV ·
   SV/fusion tables (DNA + RNA evidence) · SBS96 signatures (+ NNLS fit) · identity flags
   · oncoplot with SNV classes + CNV/fusion markers · Excel workbook · self-contained HTML report
```

## Quick start
```bash
git clone https://github.com/deppworld/mutaxplorer && cd mutaxplorer
# software: -profile conda (recommended) | docker | singularity   (envs/*.yml, docker/build_images.sh)
bash scripts/setup_resources.sh --genome hg38 --outdir resources        # one time; add --star-index for RNA

nextflow run . -profile conda --input samplesheet.csv --genome hg38 --resource_dir resources \
     --assay panel --targets panel.bed --analysis somatic --outdir results
```

### Samplesheet (CSV)
| column | meaning |
|---|---|
| `patient` | links tumour, matched normal and RNA of one patient (defaults to `sample`) |
| `sample` | unique sample ID |
| `status` | `tumor` or `normal` (default `tumor`) |
| `type` | `dna` or `rna` (default `dna`) |
| `fastq_1`, `fastq_2` | paired FASTQ — **or** `bam` (DNA only; e.g. DRAGEN / Illumina Connected Analytics) |

FASTQ and BAM rows can be mixed. A tumour with a `normal` row for the same patient is analysed tumour-normal,
otherwise tumour-only. RNA rows are used for fusion detection and linked to the patient's DNA.

## Choosing the analysis
| option | values | effect |
|---|---|---|
| `--assay` | `wgs` / `wes` / `panel` | intervals, scatter, depth/VAF thresholds, coverage bins, Mutect2 downsampling, CNVkit method, Manta mode, TMB denominator |
| `--genome` | `hg38` / `hg19` | resources, VEP assembly, MSIsensor2 models, Arriba databases |
| `--analysis` | `somatic` / `germline` / `both` | `both`: somatic on tumours, germline on normals (`--germline_samples all` for every sample) |
| `--tools` | any of `snv,cnv,sv,fusion,msi,signatures,identity` | default: all |
| `--umi` | `true` | fgbio consensus (UMI in read name, or `--umi_loc per_read --umi_len 7`) |
| `--purity` | 0–1 | rescales CNVkit integer copy numbers |
| `--pon`, `--known_sites` | VCF | Mutect2 panel of normals; BQSR |
| `--signatures` | SBS matrix | fits signature contributions (COSMIC format) |

## Tools per assay
| analysis | WGS | WES | targeted panel |
|---|---|---|---|
| read QC / trimming | fastp | fastp | fastp |
| alignment | bwa-mem2 | bwa-mem2 | bwa-mem2 |
| duplicates / UMIs | samtools markdup | samtools markdup | samtools markdup or fgbio consensus |
| base recalibration | GATK BQSR | GATK BQSR | GATK BQSR (optional) |
| coverage QC | mosdepth (genome) | mosdepth (targets) | mosdepth (targets) |
| somatic SNV/indel | GATK Mutect2 + FilterMutectCalls | same | same, downsampling off |
| germline SNV/indel | GATK HaplotypeCaller + hard filters | same | same |
| copy number | CNVkit `wgs` | CNVkit `hybrid` | CNVkit `hybrid` (`amplicon` option) |
| structural variants | Manta | Manta `--exome` | Manta `--exome` |
| DNA fusions | Manta breakpoints + gene annotation | same | same (needs intronic baits) |
| RNA fusions | STAR + Arriba | STAR + Arriba | STAR + Arriba |
| annotation | bcftools norm + Ensembl VEP | same | same |
| TMB / MSI | MutaXplorer (coding variants / Mb) / MSIsensor2 | same | same |
| signatures | SBS96 + NNLS fit | same | same (few SNVs per sample) |
| sample identity | common-SNP genotypes (bcftools) | same | same |
| aggregate QC | MultiQC | MultiQC | MultiQC |

## Outputs
```
results/
├── somatic/   vcf/ cnv/ sv/   cohort/ (all tables + Excel)   report/mutaxplorer_report.html + figures
├── germline/  vcf/ cnv/ sv/   cohort/                          report/
├── rna/fusions/               qc/ (coverage, samtools, identity)   multiqc/   pipeline_info/
```
Key tables: `variants_filtered.tsv`, `top_variants.tsv`, `top_genes.tsv`, `sample_metrics.tsv` (TMB, MSI, coverage),
`cnv_gene_calls.tsv`, `cnv_segments.tsv`, `sv_calls.tsv`, `dna_fusions.tsv`, `rna_fusions.tsv`,
`fusion_summary.tsv` (DNA / RNA / both, known-fusion flag), `sbs96_counts.tsv`, `signature_exposures.tsv`,
`sample_identity.tsv`. Every SNV keeps the first filter it failed in `variants_all_pass.tsv` (audit trail).

## Testing
```bash
bash tests/run_tests.sh      # 4 scenarios, real tools, simulated truth (SNV, indel, CNV, SV, fusion, identity)
```
`tests/benchmark_simulation.py` reproduces the recall / germline-leakage benchmark in the manuscript.

## Notes
* Imported BAMs must use the FASTA they were aligned to; MutaXplorer checks contig names, lengths and order.
* Tumour-only calling cannot fully separate rare germline variants from somatic ones; use a matched normal or PON.
* CNV from tumour-only panels relies on a pooled-normal (best) or flat reference; set `--purity` when known.
* DNA fusions from panels need baits in the fusion introns; RNA fusions need RNA-seq and a STAR index.
* TMB, MSI and CNV thresholds need local calibration before clinical use. Research use only.

## Citation and licence
Please cite the MutaXplorer manuscript (in preparation) and the tools in `CITATIONS.md`. MIT licence.
