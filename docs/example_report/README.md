# MutaXplorer example output

> **Synthetic data, not real patients.** These figures and tables come from a simulated targeted-panel
> cohort of 120 tumours on the small 18-gene test genome. They show what MutaXplorer produces, not
> biological findings.

| File | What it shows |
|---|---|
| `mutaxplorer_report.html` | The complete self-contained interactive report (download it and open it in a browser) |
| `figures/fig1_oncoplot.png` | Oncoplot: SNV/indel classes per gene and sample, with CNV and fusion markers and TMB on top |
| `figures/fig2_top_variants.png` | Most recurrent variants across the cohort |
| `figures/fig3_vaf_top_variants.png` | Variant allele fractions of the top variants |
| `figures/fig4_variant_summary.png` | Variant classes and counts per sample |
| `figures/fig5_tmb_msi.png` | TMB distribution and MSI versus TMB |
| `figures/fig6_lollipop_*.png` | Mutation positions on the protein for the most mutated genes |
| `figures/fig7_cooccurrence.png` | Co-occurrence and mutual exclusivity of mutated genes |
| `figures/cnv_genome_heatmap.png` | CNVkit segments for every sample (red gain, blue loss) |
| `figures/cnv_recurrent_genes.png` | Gene-level amplifications, losses and deep deletions |
| `figures/sv_per_sample.png` | Manta structural variants by type |
| `figures/fusions_recurrent.png` | DNA and RNA fusions with their evidence (DNA, RNA, or both) |
| `figures/signatures_sbs96.png`, `signatures_exposures.png` | SBS96 profile and NNLS signature fit |
| `tables/` | Some of the cohort tables, plus `cohort_results.xlsx` with all of them |

## Reproduce
```bash
python3 tests/simulate_test_data.py tests/data
python3 bin/gtf_to_genes.py --gtf tests/data/genes.gtf.gz --fai tests/data/genome.fa.fai --out tests/data/genes.bed
python3 tests/make_example_cohort.py --data tests/data --out example/calls
cd example && mkdir -p somatic/cohort somatic/report
D=../tests/data; B=../bin
$B/summarize_structural.py --calls-dir calls --tables somatic/cohort --report somatic/report --mode somatic \
   --samples calls/samples.tsv --gene-bed $D/genes.bed --targets $D/targets.bed --fasta $D/genome.fa \
   --known-fusions calls/known_fusions.tsv --signatures calls/example_signatures.tsv --signature-min-snvs 20 --min-dp 100
$B/aggregate_cohort.py --mode somatic --calls-dir calls --out somatic/cohort --tmb-region-mb 1.94
$B/visualize_cohort.py --mode somatic --cohort somatic/cohort --out somatic/report --top-genes 18 \
   --events somatic/cohort/oncoplot_events.tsv --extra-html somatic/report/sections.html
```
These are the same three commands the pipeline's report step runs after the real callers have finished.
