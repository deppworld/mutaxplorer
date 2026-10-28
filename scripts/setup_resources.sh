#!/usr/bin/env bash
# =============================================================================
#  MutaXplorer - download / prepare reference resources for hg38 or hg19
#
#  usage: bash scripts/setup_resources.sh --genome hg38|hg19 --outdir resources [options]
#    --fasta FILE      use your own FASTA (REQUIRED for imported BAMs: it must be the exact
#                      FASTA the BAMs were aligned with, e.g. DRAGEN/TSO500 genome.fa)
#    --targets BED     subset population resources to targets +-1 kb (WES/panel: much smaller files)
#    --bwa-index       also build the bwa-mem2 index (needs ~64 GB RAM for a human genome)
#    --star-index      also build the STAR index for RNA fusions (needs ~32 GB RAM)
#    --gencode N       GENCODE release for the gene annotation (default 46)
#    --skip-vep        do not install the VEP cache (~15-25 GB)
#
#  Produces the layout expected by:  nextflow run mutaxplorer --genome <g> --resource_dir <outdir>
#    <outdir>/<genome>/  FASTA, gnomAD, common sites, PON, (hg38: dbSNP + Mills for BQSR)
#    <outdir>/msisensor2/models_hg38 | models_hg19_GRCh37
#    <outdir>/vep_cache/
#  Needs: curl or gsutil, bcftools, samtools, tabix, git; vep_install (envs/vep.yml) for the cache.
# =============================================================================
set -euo pipefail
GENOME=""; OUT="resources"; FASTA=""; TARGETS=""; BWA=false; VEP=true; STARIDX=false; GENCODE=46
while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome) GENOME="$2"; shift 2;;  --outdir) OUT="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;    --targets) TARGETS="$2"; shift 2;;
    --bwa-index) BWA=true; shift;;    --skip-vep) VEP=false; shift;;
    --star-index) STARIDX=true; shift;; --gencode) GENCODE="$2"; shift 2;;
    -h|--help) sed -n 2,20p "$0"; exit 0;;
    *) echo "unknown option $1"; exit 1;;
  esac
done
[[ "${GENOME}" == hg38 || "${GENOME}" == hg19 ]] || { echo "--genome hg38|hg19 required"; exit 1; }
log(){ echo "[$(date '+%F %T')] $*"; }
G="${OUT}/${GENOME}"; mkdir -p "${G}" "${OUT}/download"
OUT=$(cd "${OUT}" && pwd); G="${OUT}/${GENOME}"; DL="${OUT}/download"

get(){  # get <bucket/path> <dest> : gsutil if present, else public HTTPS
  local src="$1" dst="$2"
  [[ -s "${dst}" ]] && return 0
  log "download ${src##*/}"
  if command -v gsutil >/dev/null && gsutil -q cp "gs://${src}" "${dst}" 2>/dev/null; then return 0; fi
  curl -fL --retry 3 -o "${dst}.part" "https://storage.googleapis.com/${src}" && mv "${dst}.part" "${dst}"
}
to_vcfgz(){  # compress + index (accepts .vcf or .vcf.gz)
  local src="$1" dst="$2"
  [[ -s "${dst}.tbi" ]] && return 0
  if [[ -n "${TARGETS}" ]]; then
    awk 'BEGIN{OFS="\t"} !/^(#|track|browser)/{s=$2-1000; if(s<0)s=0; print $1,s,$3+1000}' "${TARGETS}" > "${DL}/targets.pad.bed"
    bcftools view -T "${DL}/targets.pad.bed" -Oz -o "${dst}" "${src}"
  else
    bcftools view -Oz -o "${dst}" "${src}"
  fi
  tabix -f -p vcf "${dst}"
}

# ---- 1. FASTA --------------------------------------------------------------
if [[ "${GENOME}" == hg38 ]]; then REF="${G}/Homo_sapiens_assembly38.fasta"; else REF="${G}/hg19.fa"; fi
if [[ -n "${FASTA}" ]]; then
  ln -sf "$(cd "$(dirname "${FASTA}")" && pwd)/$(basename "${FASTA}")" "${REF}"
elif [[ ! -s "${REF}" ]]; then
  if [[ "${GENOME}" == hg38 ]]; then
    get genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta "${REF}"
  else
    log "download UCSC hg19"; curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gunzip > "${REF}"
  fi
fi
[[ -s "${REF}.fai" ]] || samtools faidx "${REF}"
[[ -s "${REF%.*}.dict" ]] || samtools dict "${REF}" -o "${REF%.*}.dict"
CHR=""; head -1 "${REF}.fai" | grep -q '^chr' && CHR="chr"

# ---- 2. Mutect2 / BQSR resources ------------------------------------------
if [[ "${GENOME}" == hg38 ]]; then
  B=gatk-best-practices/somatic-hg38
  for f in af-only-gnomad.hg38.vcf.gz small_exac_common_3.hg38.vcf.gz 1000g_pon.hg38.vcf.gz; do
    get "${B}/${f}" "${DL}/${f}" && get "${B}/${f}.tbi" "${DL}/${f}.tbi" || log "WARNING could not fetch ${f}"
    [[ -s "${DL}/${f}" ]] && to_vcfgz "${DL}/${f}" "${G}/${f}"
  done
  R=genomics-public-data/references/hg38/v0
  get "${R}/Homo_sapiens_assembly38.dbsnp138.vcf" "${DL}/dbsnp138.vcf" && to_vcfgz "${DL}/dbsnp138.vcf" "${G}/Homo_sapiens_assembly38.dbsnp138.vcf.gz" || log "WARNING dbSNP not fetched (BQSR)"
  get "${R}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" "${DL}/mills.vcf.gz" && to_vcfgz "${DL}/mills.vcf.gz" "${G}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" || log "WARNING Mills not fetched (BQSR)"
else
  # hg19: GATK publishes these for b37 (no 'chr'); same coordinates for 1-22,X,Y -> rename, drop MT
  B=gatk-best-practices/somatic-b37
  for c in $(seq 1 22) X Y; do echo -e "${c}\t${CHR}${c}"; done > "${DL}/b37_rename.txt"
  KEEP_B37=$(cut -f1 "${DL}/b37_rename.txt" | paste -sd,)
  b37(){  # b37 <name> <hg19 output>
    local n="$1" o="$2"; [[ -s "${o}.tbi" ]] && return 0
    get "${B}/${n}" "${DL}/${n}" || { log "WARNING could not fetch ${n} - download it manually into ${DL}/"; return 0; }
    bcftools view -t "${KEEP_B37}" -Ou "${DL}/${n}" 2>/dev/null \
      | bcftools annotate --rename-chrs "${DL}/b37_rename.txt" -Ob -o "${DL}/${n}.renamed.bcf"
    bcftools reheader -f "${REF}.fai" -o "${DL}/${n}.rh.bcf" "${DL}/${n}.renamed.bcf"
    bcftools sort -Oz -o "${DL}/${n}.hg19.vcf.gz" "${DL}/${n}.rh.bcf" && tabix -f -p vcf "${DL}/${n}.hg19.vcf.gz"
    to_vcfgz "${DL}/${n}.hg19.vcf.gz" "${o}"; rm -f "${DL}/${n}".*.bcf
  }
  b37 af-only-gnomad.raw.sites.vcf "${G}/af-only-gnomad.hg19.vcf.gz"
  b37 small_exac_common_3.vcf       "${G}/small_exac_common_3.hg19.vcf.gz"
  b37 Mutect2-exome-panel.vcf       "${G}/Mutect2-exome-panel.hg19.vcf.gz"
  log "hg19 BQSR known sites are not downloaded automatically: pass --known_sites dbsnp.vcf.gz,mills.vcf.gz (hg19, indexed)"
fi

# ---- 3. gene annotation (GENCODE) + known fusions (from Arriba's database) ---------
GTF="${G}/gencode.annotation.gtf.gz"
if [[ ! -s "${GTF}" ]]; then
  if [[ "${GENOME}" == hg38 ]]; then U="release_${GENCODE}/gencode.v${GENCODE}.annotation.gtf.gz"
  else U="release_${GENCODE}/GRCh37_mapping/gencode.v${GENCODE}lift37.annotation.gtf.gz"; fi
  log "download GENCODE ${GENCODE} annotation"
  curl -fL --retry 3 -o "${GTF}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/${U}" \
    || { rm -f "${GTF}"; log "WARNING GENCODE download failed - set --gtf manually"; }
fi
if [[ ! -s "${G}/known_fusions.tsv.gz" ]] && command -v arriba >/dev/null; then
  B=$([[ "${GENOME}" == hg38 ]] && echo hg38_GRCh38 || echo hg19_hs37d5_GRCh37)
  KF=$(ls "$(dirname "$(command -v arriba)")"/../var/lib/arriba/known_fusions_${B}_*.tsv.gz 2>/dev/null | tail -1)
  [[ -n "${KF}" ]] && cp "${KF}" "${G}/known_fusions.tsv.gz" && log "known fusions: $(basename "${KF}")"
fi
[[ -s "${G}/known_fusions.tsv.gz" ]] || log "NOTE known-fusion list not installed (activate envs/arriba.yml and re-run, or pass --known_fusions)"

# ---- 3b. MSIsensor2 models ----------------------------------------------------
if [[ ! -d "${OUT}/msisensor2/models_hg38" ]]; then
  log "MSIsensor2 models"; rm -rf "${OUT}/msisensor2"
  git clone -q --depth 1 https://github.com/niu-lab/msisensor2.git "${OUT}/msisensor2"
fi

# ---- 4. VEP cache --------------------------------------------------------------
if ${VEP}; then
  ASM=$([[ "${GENOME}" == hg38 ]] && echo GRCh38 || echo GRCh37)
  if compgen -G "${OUT}/vep_cache/homo_sapiens*/*_${ASM}" >/dev/null; then log "VEP ${ASM} cache present"
  elif command -v vep_install >/dev/null; then
    log "VEP ${ASM} cache (~15-25 GB)"; mkdir -p "${OUT}/vep_cache"
    vep_install -a c -s homo_sapiens -y "${ASM}" -c "${OUT}/vep_cache" --CONVERT --NO_UPDATE
  else
    log "WARNING vep_install not found: conda env create -f envs/vep.yml && conda activate mutaxplorer-vep, then re-run"
  fi
fi

# ---- 5. bwa-mem2 index ---------------------------------------------------------
if ${BWA} && [[ ! -s "${G}/bwamem2/$(basename "${REF}").bwt.2bit.64" ]]; then
  log "bwa-mem2 index (long, ~64 GB RAM)"; mkdir -p "${G}/bwamem2"
  ln -sf "${REF}" "${G}/bwamem2/$(basename "${REF}")"; bwa-mem2 index "${G}/bwamem2/$(basename "${REF}")"
fi
# ---- 6. STAR index for RNA fusions -------------------------------------------------
if ${STARIDX} && [[ ! -s "${G}/star/SA" ]] && [[ -s "${GTF}" ]]; then
  log "STAR index (long, ~32 GB RAM)"; mkdir -p "${G}/star"
  zcat "${GTF}" > "${G}/star/genes.gtf"
  STAR --runMode genomeGenerate --runThreadN 8 --genomeDir "${G}/star" --genomeFastaFiles "${REF}" \
       --sjdbGTFfile "${G}/star/genes.gtf" --sjdbOverhang 149
fi
log "Mutational signatures: download the COSMIC SBS matrix for ${GENOME} (cancer.sanger.ac.uk/signatures/downloads) and pass --signatures"
log "Done. Use:  --genome ${GENOME} --resource_dir ${OUT}$(${BWA} && echo " --bwa_index ${G}/bwamem2")$(${STARIDX} && echo " --star_index ${G}/star")"
