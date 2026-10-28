# Installing MutaXplorer

MutaXplorer is a Nextflow pipeline. You install **Nextflow** plus **one way of providing the
bioinformatics tools** (conda/mamba, Docker or Singularity). Nextflow then sets up every tool
automatically the first time the pipeline runs. You do not install GATK, CNVkit, Manta and the rest by hand.

## 1. System requirements

| | Minimum |
|---|---|
| Operating system | Linux (x86-64) or macOS. **Windows: use WSL2** (Ubuntu) — the tools do not run natively on Windows |
| CPU / memory, targeted panel | 8 cores, 32 GB RAM |
| CPU / memory, exome | 16 cores, 64 GB RAM |
| CPU / memory, genome | 32 cores, 128 GB RAM |
| Building indexes (one time) | bwa-mem2 index ~64 GB RAM; STAR index (RNA fusions) ~32 GB RAM |
| Disk for reference resources | ~100 GB (hg38 or hg19 FASTA + indexes, VEP cache ~20 GB, gnomAD resources) |
| Internet | needed for the first run (downloads software and resources) |

Windows users: open PowerShell as administrator, run `wsl --install`, restart, open "Ubuntu",
and do everything below inside Ubuntu.

## 2. Install the prerequisites

### Java and Nextflow
Nextflow needs Java 17 or newer.
```bash
# Java (Ubuntu/WSL)
sudo apt update && sudo apt install -y openjdk-17-jre-headless git curl
# Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version          # should print 24.04 or newer
```

### Choose ONE way to provide the tools

**A. Conda / Mamba (simplest; no administrator rights needed)**
```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b && ~/miniforge3/bin/conda init && exec bash
```
Then run the pipeline with `-profile conda` (or `-profile mamba`, faster).

**B. Docker**: install Docker Desktop (Windows/macOS) or Docker Engine (Linux), then build the images once:
```bash
bash docker/build_images.sh
```
Run with `-profile docker`.

**C. Singularity / Apptainer** (HPC clusters): ask your HPC team to load it, build the images as in B,
and run with `-profile singularity` or `-profile apptainer`.

## 3. Get MutaXplorer
```bash
git clone https://github.com/deppworld/mutaxplorer.git
cd mutaxplorer
```

## 4. Check the installation (built-in test, ~10 minutes)
```bash
conda env create -f envs/mutaxplorer.yml && conda activate mutaxplorer   # for the test data generator
python3 tests/simulate_test_data.py tests/data
nextflow run . -profile test,conda          # or test,docker / test,singularity
```
The run should end with `Succeeded` and write `results/somatic/report/mutaxplorer_report.html`.
The test uses a small simulated genome and a mock of VEP, so no large downloads are needed.

## 5. Download the reference resources (one time)
```bash
conda env create -f envs/vep.yml          # VEP, needed to install its cache
conda env create -f envs/arriba.yml       # provides the known-fusion list
conda activate mutaxplorer-vep
bash scripts/setup_resources.sh --genome hg38 --outdir resources
# options: --genome hg19 | --fasta <your.fa> (REQUIRED for DRAGEN/vendor BAMs: the FASTA they were aligned to)
#          --targets panel.bed (smaller resource files for WES/panels) | --bwa-index | --star-index
```
This downloads the reference genome, GATK germline/contamination resources, GENCODE gene annotation,
MSIsensor2 models and the Ensembl VEP cache. For mutational-signature fitting, download the COSMIC SBS
matrix for your genome build from cancer.sanger.ac.uk/signatures/downloads and pass `--signatures`.

## 6. Run on your data
```bash
nextflow run . -profile conda \
    --input samplesheet.csv --genome hg38 --resource_dir resources \
    --assay panel --targets panel.bed --analysis somatic --outdir results
```
See `README.md` for the samplesheet format and all options, and `assets/` for examples.
On a SLURM cluster add `-profile conda,slurm`. Add `-resume` to continue an interrupted run.

## Tools installed automatically

| Tool | Version | Used for | Environment file |
|---|---|---|---|
| Nextflow | ≥ 24.04 | workflow engine | install yourself (step 2) |
| fastp | ≥ 0.23 | read QC and trimming | `envs/mutaxplorer.yml` |
| bwa-mem2 | ≥ 2.2.1 | alignment | `envs/mutaxplorer.yml` |
| samtools, bcftools, htslib | ≥ 1.20 | BAM/VCF handling, normalisation, identity genotypes | `envs/mutaxplorer.yml` |
| GATK | ≥ 4.5 | Mutect2, HaplotypeCaller, BQSR, filtering | `envs/mutaxplorer.yml` |
| fgbio | ≥ 2.2 | UMI consensus | `envs/mutaxplorer.yml` |
| mosdepth | ≥ 0.3.6 | coverage | `envs/mutaxplorer.yml` |
| MSIsensor2 | latest | MSI | `envs/mutaxplorer.yml` |
| MultiQC | ≥ 1.20 | QC report | `envs/mutaxplorer.yml` |
| Python 3 + pandas, numpy, scipy, matplotlib, openpyxl | ≥ 3.10 | cohort tables and report | `envs/mutaxplorer.yml` |
| CNVkit | ≥ 0.9.10 | copy number | `envs/cnvkit.yml` |
| Manta | ≥ 1.6.0 | structural variants, DNA fusions | `envs/manta.yml` |
| STAR, Arriba | ≥ 2.7.10, ≥ 2.4 | RNA fusions | `envs/arriba.yml` |
| Ensembl VEP | ≥ 110 | variant annotation | `envs/vep.yml` |

With `-profile conda`, Nextflow creates these four environments automatically the first time (10–30 minutes).

## Troubleshooting
- **`nextflow: command not found`**: Nextflow is not on your PATH (step 2).
- **Conda environment creation is slow or fails**: use `-profile mamba`, or create the environments yourself
  with `conda env create -f envs/<file>.yml`.
- **"contigs do not match the reference" for BAM input**: use `--fasta` with the exact FASTA the BAMs
  were aligned to (for DRAGEN TSO500, the hg19 `genome.fa` from the TSO500 resource bundle).
- **Out of memory**: lower the limits with `--max_cpus 8 --max_memory 32.GB`, or run on a larger machine.
