#!/usr/bin/env nextflow
/*
 * ============================================================================
 *  MutaXplorer - comprehensive DNA + RNA variant exploration
 *  WGS, WES or targeted panels | FASTQ or BAM | hg19 or hg38
 *  SNV/indel (somatic + germline), copy number, structural variants,
 *  DNA + RNA fusions, TMB, MSI, mutational signatures, sample identity,
 *  QC and an interactive cohort report.
 * ============================================================================
 */

include { PREPARE_REFERENCE; BWA_INDEX; MAKE_INTERVALS; GENE_BED } from './modules/local/reference'
include { FASTP; ALIGN; MARK_DUPLICATES; UMI_CONSENSUS; BQSR; PREPARE_BAM } from './modules/local/preprocess'
include { MOSDEPTH; SAMTOOLS_STATS; MULTIQC; IDENTITY            } from './modules/local/qc'
include { MUTECT2; FILTER_MUTECT; MSISENSOR2                     } from './modules/local/somatic'
include { HAPLOTYPECALLER; GERMLINE_FILTER                       } from './modules/local/germline'
include { NORMALIZE as NORMALIZE_SOM; NORMALIZE as NORMALIZE_GERM } from './modules/local/annotate'
include { VEP as VEP_SOM; VEP as VEP_GERM                        } from './modules/local/annotate'
include { CNVKIT_REFERENCE as CNV_REF_SOM; CNVKIT_REFERENCE as CNV_REF_GERM } from './modules/local/cnv'
include { CNVKIT_CALL as CNV_CALL_SOM; CNVKIT_CALL as CNV_CALL_GERM } from './modules/local/cnv'
include { MANTA as MANTA_SOM; MANTA as MANTA_GERM                 } from './modules/local/sv'
include { STAR_INDEX; ARRIBA                                     } from './modules/local/rna'
include { COHORT_REPORT as REPORT_SOM; COHORT_REPORT as REPORT_GERM } from './modules/local/report'

// ---------------------------------------------------------------- helpers
// explicit --param wins; otherwise the genome's file under --resource_dir, if it exists
def res(key) {
    if (params[key]) return params[key]
    def rel = params.genomes[params.genome]?.get(key)
    if (!params.resource_dir || !rel) return null
    def paths = rel.tokenize(',').collect { "${params.resource_dir}/${it}".toString() }
    return paths.every { file(it).exists() } ? paths.join(',') : null
}

def optionalFile(p) { p ? file(p, checkIfExists: true) : [] }

def optionalVcf(p) {
    if (!p) return [[], []]
    def v = file(p, checkIfExists: true)
    def idx = file("${p}.tbi").exists() ? file("${p}.tbi") : file("${p}.idx")
    if (!idx.exists()) error "Index (.tbi/.idx) missing for ${p}"
    return [v, idx]
}

def tools() { params.tools.toString().toLowerCase().tokenize(',')*.trim() }
def useTool(t) { t in tools() }

def validateParams() {
    def errors = []
    def known = ['snv', 'cnv', 'sv', 'fusion', 'msi', 'signatures', 'identity']
    if (!params.input)  errors << "--input <samplesheet.csv> is required"
    if (!res('fasta'))  errors << "--fasta (or --genome with --resource_dir) is required"
    if (!(params.genome in ['hg19', 'hg38']))                  errors << "--genome must be hg19 or hg38"
    if (!(params.assay in ['wgs', 'wes', 'panel']))            errors << "--assay must be wgs, wes or panel"
    if (!(params.analysis in ['somatic', 'germline', 'both'])) errors << "--analysis must be somatic, germline or both"
    if (params.assay != 'wgs' && !params.targets)              errors << "--targets BED is required for --assay ${params.assay}"
    def bad = tools() - known
    if (bad) errors << "unknown --tools ${bad}; choose from ${known.join(',')}"
    if (useTool('fusion') && !res('gtf')) log.warn "--gtf not set: DNA fusion annotation and RNA fusion calling are skipped"
    if (useTool('cnv') && !res('gtf'))    log.warn "--gtf not set: CNV segments are reported without gene-level calls"
    if (errors) error "Parameter error(s):\n  - " + errors.join("\n  - ")
}

def parseSamplesheet(path) {
    channel.fromPath(path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample) error "Samplesheet row without 'sample': ${row}"
            def status = (row.status ?: 'tumor').toLowerCase()
            def type   = (row.type ?: 'dna').toLowerCase()
            if (!(status in ['tumor', 'normal'])) error "status must be tumor or normal (sample ${row.sample})"
            if (!(type in ['dna', 'rna']))        error "type must be dna or rna (sample ${row.sample})"
            def meta = [id: row.sample, patient: row.patient ?: row.sample, status: status, type: type, lane: row.lane ?: '1']
            if (row.bam) {
                if (type == 'rna') error "Sample ${row.sample}: RNA input must be FASTQ"
                return [meta + [input: 'bam'], file(row.bam, checkIfExists: true), [], []]
            }
            if (row.fastq_1 && row.fastq_2) {
                return [meta + [input: 'fastq'], [], file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
            }
            error "Sample ${row.sample}: give either 'bam' or 'fastq_1' + 'fastq_2'"
        }
}

// ---------------------------------------------------------------- workflow
workflow {
    validateParams()

    log.info """
    ==================================================
     M U T A X P L O R E R   v${workflow.manifest.version}
    ==================================================
     input     : ${params.input}
     genome    : ${params.genome}   fasta: ${res('fasta')}
     assay     : ${params.assay}    targets: ${params.targets ?: '-'}
     analysis  : ${params.analysis}
     tools     : ${tools().join(', ')}
     UMI       : ${params.umi}      BQSR: ${res('known_sites') ? 'yes' : 'no'}
     annotation: ${res('vep_cache') ? 'VEP ' + params.vep_assembly : 'none (no --vep_cache)'}   genes: ${res('gtf') ?: '-'}
     outdir    : ${params.outdir}
    ==================================================
    """.stripIndent()

    // ---- reference, intervals, genes --------------------------------------------
    PREPARE_REFERENCE(file(res('fasta'), checkIfExists: true))
    ref = PREPARE_REFERENCE.out.ref.collect()
    fai = ref.map { it[1] }
    MAKE_INTERVALS(fai, optionalFile(params.targets), optionalFile(params.exclude_bed), params.scatter_count)
    intervals_all = MAKE_INTERVALS.out.all.collect()
    chunks        = MAKE_INTERVALS.out.chunks.flatten()
    gtf = res('gtf') ? file(res('gtf'), checkIfExists: true) : null
    genes = gtf ? GENE_BED(gtf, fai).bed.collect() : channel.value([])

    samples = parseSamplesheet(params.input)
    dna = samples.filter { it[0].type == 'dna' }
    rna = samples.filter { it[0].type == 'rna' }.map { m, b, r1, r2 -> [m, r1, r2] }
    fq_in  = dna.filter { it[0].input == 'fastq' }.map { m, b, r1, r2 -> [m, r1, r2] }
    bam_in = dna.filter { it[0].input == 'bam' }.map { m, b, r1, r2 -> [m, b] }

    // sample table for the reports (id, patient, status, type)
    samples_tsv = samples.map { it[0] }.map { m -> "${m.id}\t${m.patient}\t${m.status}\t${m.type}" }
        .collectFile(name: 'samples.tsv', newLine: true, sort: true).collect()

    // ---- FASTQ -> analysis-ready BAM ------------------------------------------------
    bwa_index = params.bwa_index ? channel.value(file(params.bwa_index, checkIfExists: true))
                                 : BWA_INDEX(fq_in.first().map { file(res('fasta')) }).index.collect()
    FASTP(fq_in)
    ALIGN(FASTP.out.reads, ref, bwa_index)
    if (params.umi) {
        UMI_CONSENSUS(ALIGN.out.bam, ref, bwa_index)
        fq_bams = UMI_CONSENSUS.out.bam
        dedup_mqc = UMI_CONSENSUS.out.metrics
    } else {
        MARK_DUPLICATES(ALIGN.out.bam)
        fq_bams = MARK_DUPLICATES.out.bam
        dedup_mqc = MARK_DUPLICATES.out.metrics
    }
    if (res('known_sites')) {
        known = channel.fromPath(res('known_sites').tokenize(',').collectMany { [it, "${it}.tbi"] }, checkIfExists: true).collect()
        BQSR(fq_bams, ref, known, intervals_all)
        fq_bams = BQSR.out.bam
    }
    PREPARE_BAM(bam_in, ref)
    bams = fq_bams.mix(PREPARE_BAM.out.bam)        // [meta, bam, bai]

    // ---- QC + sample identity ------------------------------------------------------------
    MOSDEPTH(bams, intervals_all, params.coverage_thresholds)
    SAMTOOLS_STATS(bams)
    identity = channel.empty()
    if (useTool('identity') && res('common_sites')) {
        IDENTITY(bams.map { it[1] }.collect(), bams.map { it[2] }.collect(), ref, channel.value(optionalVcf(res('common_sites'))))
        identity = IDENTITY.out.gtcheck
    }

    known_fusions = res('known_fusions') ? file(res('known_fusions'), checkIfExists: true) : []
    signatures    = res('signatures') ? file(res('signatures'), checkIfExists: true) : []

    // ---- SOMATIC -----------------------------------------------------------------------------
    if (params.analysis in ['somatic', 'both']) {
        normals_by_patient = bams.filter { it[0].status == 'normal' }
            .map { m, b, i -> [m.patient, [m, b, i]] }
            .toList()
            .map { lst -> lst.collectEntries { [(it[0]): it[1]] } }
        pairs = bams.filter { it[0].status == 'tumor' }
            .combine(normals_by_patient)
            .map { m, b, i, nmap ->
                def n = nmap[m.patient]
                n ? [m + [mode: 'somatic', pair: 'tumor_normal', normal_id: n[0].id], b, i, n[1], n[2]]
                  : [m + [mode: 'somatic', pair: 'tumor_only'], b, i, [], []]
            }
        tumors = pairs.map { m, b, i, nb, ni -> [m, b, i] }
        som_files = channel.empty()

        // SNV / indel
        MUTECT2(pairs.combine(chunks), ref,
                channel.value(optionalVcf(res('germline_resource'))), channel.value(optionalVcf(res('pon'))))
        FILTER_MUTECT(MUTECT2.out.calls.groupTuple().join(pairs), ref,
                      channel.value(optionalVcf(res('common_sites'))), intervals_all)
        NORMALIZE_SOM(FILTER_MUTECT.out.vcf, ref)
        som_files = som_files.mix(NORMALIZE_SOM.out.vcf.map { m, v, t -> v })
        if (res('vep_cache')) {
            VEP_SOM(NORMALIZE_SOM.out.vcf, ref, file(res('vep_cache'), checkIfExists: true), params.vep_assembly)
            som_files = som_files.mix(VEP_SOM.out.vcf.map { m, v -> v })
        }
        // MSI
        if (useTool('msi') && res('msisensor2_models')) {
            MSISENSOR2(tumors, file(res('msisensor2_models'), checkIfExists: true))
            som_files = som_files.mix(MSISENSOR2.out.msi.map { m, f -> f })
        }
        // copy number: reference from the run's normals (flat if none)
        if (useTool('cnv')) {
            normal_bams = bams.filter { it[0].status == 'normal' }
            CNV_REF_SOM(normal_bams.map { it[1] }.collect().ifEmpty([]), normal_bams.map { it[2] }.collect().ifEmpty([]),
                        ref, optionalFile(params.targets) ? MAKE_INTERVALS.out.all.collect() : [])
            CNV_CALL_SOM(tumors, CNV_REF_SOM.out.cnn.collect())
            som_files = som_files.mix(CNV_CALL_SOM.out.calls.flatMap { m, c, r -> [c, r] })
        }
        // structural variants (DNA fusions are derived from them in the report)
        if (useTool('sv') || useTool('fusion')) {
            MANTA_SOM(pairs, ref, intervals_all, 'somatic')
            som_files = som_files.mix(MANTA_SOM.out.vcf.map { m, v, t -> v })
        }
        // RNA fusions
        if (useTool('fusion') && gtf) {
            star = params.star_index ? channel.value(file(params.star_index, checkIfExists: true))
                                     : STAR_INDEX(rna.first().map { 1 }.combine(ref).map { it[1..3] }, gtf).index.collect()
            ARRIBA(rna, ref, gtf, star)
            som_files = som_files.mix(ARRIBA.out.fusions.map { m, f -> f })
        }
        tumor_ids = tumors.map { it[0].id }.collect()
        som_cov = MOSDEPTH.out.cov.combine(tumor_ids.map { [it] })
            .filter { m, s, t, ids -> m.id in ids }.flatMap { m, s, t, ids -> [s, t] }
        REPORT_SOM('somatic', som_files.mix(som_cov, identity).collect(), samples_tsv, ref, genes,
                   intervals_all, known_fusions, useTool('signatures') ? signatures : [])
    }

    // ---- GERMLINE --------------------------------------------------------------------------------
    if (params.analysis in ['germline', 'both']) {
        germ_bams = bams.filter { params.analysis == 'germline' || params.germline_samples == 'all' || it[0].status == 'normal' }
            .map { m, b, i -> [m + [mode: 'germline'], b, i] }
        germ_files = channel.empty()
        HAPLOTYPECALLER(germ_bams.combine(chunks), ref)
        GERMLINE_FILTER(HAPLOTYPECALLER.out.calls.groupTuple(), ref)
        NORMALIZE_GERM(GERMLINE_FILTER.out.vcf, ref)
        germ_files = germ_files.mix(NORMALIZE_GERM.out.vcf.map { m, v, t -> v })
        if (res('vep_cache')) {
            VEP_GERM(NORMALIZE_GERM.out.vcf, ref, file(res('vep_cache'), checkIfExists: true), params.vep_assembly)
            germ_files = germ_files.mix(VEP_GERM.out.vcf.map { m, v -> v })
        }
        if (useTool('cnv')) {   // germline CNV against a flat reference
            CNV_REF_GERM([], [], ref, optionalFile(params.targets) ? MAKE_INTERVALS.out.all.collect() : [])
            CNV_CALL_GERM(germ_bams, CNV_REF_GERM.out.cnn.collect())
            germ_files = germ_files.mix(CNV_CALL_GERM.out.calls.flatMap { m, c, r -> [c, r] })
        }
        if (useTool('sv')) {
            MANTA_GERM(germ_bams.map { m, b, i -> [m, b, i, [], []] }, ref, intervals_all, 'germline')
            germ_files = germ_files.mix(MANTA_GERM.out.vcf.map { m, v, t -> v })
        }
        germ_ids = germ_bams.map { it[0].id }.collect()
        germ_cov = MOSDEPTH.out.cov.combine(germ_ids.map { [it] })
            .filter { m, s, t, ids -> m.id in ids }.flatMap { m, s, t, ids -> [s, t] }
        REPORT_GERM('germline', germ_files.mix(germ_cov, identity).collect(), samples_tsv, ref, genes,
                    intervals_all, [], [])
    }

    // ---- MultiQC -----------------------------------------------------------------------------------
    MULTIQC(FASTP.out.json.mix(dedup_mqc, MOSDEPTH.out.mqc.flatten(), SAMTOOLS_STATS.out.mqc.flatten()).collect())
}
