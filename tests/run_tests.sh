#!/usr/bin/env bash
# =============================================================================
#  MutaXplorer test suite - real tools on simulated data with known truth (~15 min, 2 CPUs)
#   1 panel | somatic + germline | FASTQ tumour/normal + RNA + imported BAM (no read group)
#           checks SNV/indel, germline, CNV (amp / deep del / loss), EML4--ALK (DNA + RNA), identity
#   2 panel | somatic            | UMI FASTQ (fgbio consensus)
#   3 WGS   | germline           | FASTQ, bwa-mem2 index built by the pipeline, germline SV (2-kb DEL)
#   4 WES   | somatic + germline | hg19 via --resource_dir layout, BQSR
#  Annotation uses tests/bin/vep (a MOCK of VEP) so the report path runs without the 20 GB cache.
#  usage: bash tests/run_tests.sh [workdir]
#  needs on PATH: envs/mutaxplorer.yml tools + wgsim, and cnvkit.py, configManta.py, STAR + arriba
#  (or run each scenario with -profile conda).
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
W="$(mkdir -p "${1:-${HERE}/test_runs}" && cd "${1:-${HERE}/test_runs}" && pwd)"
D="${W}/data"
export PATH="${HERE}/tests/bin:${PATH}" NXF_ANSI_LOG=false MUTAXPLORER_MOCK_GNOMAD="${D}/gnomad.vcf.gz"
fail=0
[[ -s "${D}/truth.tsv" ]] || python3 "${HERE}/tests/simulate_test_data.py" "${D}" || exit 1

RD="${W}/resdir"; mkdir -p "${RD}/hg19"
ln -sf "${D}/genome.fa" "${RD}/hg19/hg19.fa"; ln -sf "${D}/genes.gtf.gz" "${RD}/hg19/gencode.annotation.gtf.gz"
for f in gnomad:af-only-gnomad common:small_exac_common_3; do
  ln -sf "${D}/${f%%:*}.vcf.gz" "${RD}/hg19/${f#*:}.hg19.vcf.gz"; ln -sf "${D}/${f%%:*}.vcf.gz.tbi" "${RD}/hg19/${f#*:}.hg19.vcf.gz.tbi"
done

run(){
  local name="$1"; shift
  mkdir -p "${W}/${name}"
  ( cd "${W}/${name}" && nextflow run "${HERE}" -c "${HERE}/tests/ci_limits.config" --outdir results -resume "$@" > run.log 2>&1 )
  local rc=$?
  echo "== ${name}: nextflow exit ${rc}"; [[ ${rc} -eq 0 ]] || { grep -A20 "Caused by" "${W}/${name}/run.log" | head -30; fail=1; }
}
check(){ python3 "${HERE}/tests/check_results.py" "${D}/truth.tsv" "${W}/$1/results" "${@:2}" || fail=1; }

CNV="P1_T:ERBB2:Amplification,P1_T:CDKN2A:Deep deletion,P1_T:PTEN:Loss,P2_T:ERBB2:Amplification,P2_T:CDKN2A:Deep deletion,P2_T:PTEN:Loss"
COMMON=(--fasta "${D}/genome.fa" --targets "${D}/targets.bed" --gtf "${D}/genes.gtf.gz" --vep_cache "${D}/vep_cache"
        --germline_resource "${D}/gnomad.vcf.gz" --common_sites "${D}/common.vcf.gz" --arriba_blacklist none)

run panel  --input "${D}/samplesheet.csv" "${COMMON[@]}" --bwa_index "${D}/bwamem2" --assay panel --analysis both \
           --known_fusions "${D}/known_fusions.tsv" --signatures "${D}/toy_signatures.tsv" --signature_min_snvs 3
check panel --somatic P1_T,P2_T --tumor-only P2_T --germline P1_N --cnv "${CNV}" \
            --fusion "P1:EML4--ALK:DNA + RNA,P2:EML4--ALK:DNA" --identity "P1_T:P1_N:OK,P1_N:P2_T:UNEXPECTED"
run umi    --input "${D}/samplesheet_umi.csv" "${COMMON[@]}" --bwa_index "${D}/bwamem2" --assay panel --analysis somatic \
           --umi true --tools snv,cnv,sv
check umi  --somatic P3_T --tumor-only P3_T --cnv "P3_T:ERBB2:Amplification,P3_T:CDKN2A:Deep deletion"
run wgs    --input "${D}/samplesheet_wgs.csv" --fasta "${D}/genome.fa" --gtf "${D}/genes.gtf.gz" --vep_cache "${D}/vep_cache" \
           --assay wgs --analysis germline --scatter_count 4 --tools snv,cnv,sv
check wgs  --germline WGS_N --germline-sv "WGS_N:DEL:chr3:115000"
run wes    --input "${D}/samplesheet.csv" --genome hg19 --resource_dir "${RD}" --targets "${D}/targets.bed" \
           --vep_cache "${D}/vep_cache" --bwa_index "${D}/bwamem2" --assay wes --analysis both \
           --known_sites "${D}/known_sites.vcf.gz" --arriba_blacklist none
check wes  --somatic P1_T,P2_T --tumor-only P2_T --germline P1_N --cnv "${CNV}" --fusion "P1:EML4--ALK:DNA + RNA"

for r in panel umi wgs wes; do
  for m in somatic germline; do
    [[ -d "${W}/${r}/results/${m}" ]] && { [[ -s "${W}/${r}/results/${m}/report/mutaxplorer_report.html" ]] \
      && echo "  [PASS] ${r}: ${m} report written" || { echo "  [FAIL] ${r}: ${m} report missing"; fail=1; }; }
  done
done
[[ ${fail} -eq 0 ]] && echo "MUTAXPLORER TEST SUITE PASSED" || { echo "MUTAXPLORER TEST SUITE FAILED"; exit 1; }
