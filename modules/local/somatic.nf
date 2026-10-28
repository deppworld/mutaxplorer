// MutaXplorer - somatic SNV/indel calling (GATK4 Mutect2: tumor-only or tumor-normal) + MSI

process MUTECT2 {
    tag "${meta.id}:${interval.baseName.tokenize('.').last()}"
    label 'process_medium'

    input:
    tuple val(meta), path(tbam), path(tbai), path(nbam), path(nbai), path(interval)
    tuple path(fasta), path(fai), path(dict)
    tuple path(gnomad), path(gnomad_tbi)
    tuple path(pon), path(pon_tbi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), path("${prefix}.vcf.gz.stats"), path("${prefix}.f1r2.tar.gz"), emit: calls

    script:
    prefix = "${meta.id}.${interval.baseName.tokenize('.').last()}"
    def normal  = nbam ? "-I ${nbam} -normal \$(samtools view -H ${nbam} | grep -m1 '^@RG' | tr '\\t' '\\n' | awk -F: '\$1==\"SM\"{print \$2}')" : ''
    def germ    = gnomad ? "--germline-resource ${gnomad}" : ''
    def panel   = pon ? "--panel-of-normals ${pon}" : ''
    """
    gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" Mutect2 \\
        -R ${fasta} -I ${tbam} ${normal} -L ${interval} \\
        ${germ} ${panel} \\
        --max-reads-per-alignment-start ${params.max_reads_per_start} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz -O ${prefix}.vcf.gz ${params.mutect2_args}
    """
}

process FILTER_MUTECT {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcfs), path(tbis), path(stats), path(f1r2), path(tbam), path(tbai), path(nbam), path(nbai)
    tuple path(fasta), path(fai), path(dict)
    tuple path(common), path(common_tbi)
    path intervals

    output:
    tuple val(meta), path("${meta.id}.mutect2.filtered.vcf.gz"), path("${meta.id}.mutect2.filtered.vcf.gz.tbi"), emit: vcf
    path "${meta.id}.contamination.table", optional: true, emit: contamination

    script:
    def mem = task.memory.toGiga() - 1
    def ins = (vcfs instanceof List ? vcfs : [vcfs]).collect { "-I ${it}" }.join(' ')
    def sts = (stats instanceof List ? stats : [stats]).collect { "-stats ${it}" }.join(' ')
    def f1s = (f1r2 instanceof List ? f1r2 : [f1r2]).collect { "-I ${it}" }.join(' ')
    """
    G() { gatk --java-options "-Xmx${mem}g" "\$@"; }
    G MergeVcfs ${ins} -O merged.vcf.gz
    G MergeMutectStats ${sts} -O merged.vcf.gz.stats
    G LearnReadOrientationModel ${f1s} -O read-orientation-model.tar.gz
    EXTRA=""
    if [[ -n "${common}" ]]; then
        if G GetPileupSummaries -I ${tbam} -V ${common} -L ${common} -L ${intervals} \\
               --interval-set-rule INTERSECTION -O tumor.pileups.table; then
            MATCHED=""
            if [[ -n "${nbam}" ]] && G GetPileupSummaries -I ${nbam} -V ${common} -L ${common} -L ${intervals} \\
                   --interval-set-rule INTERSECTION -O normal.pileups.table; then
                MATCHED="-matched normal.pileups.table"
            fi
            if G CalculateContamination -I tumor.pileups.table \$MATCHED \\
                   -O ${meta.id}.contamination.table --tumor-segmentation segments.table; then
                EXTRA="--contamination-table ${meta.id}.contamination.table --tumor-segmentation segments.table"
            fi
        fi
        [[ -n "\$EXTRA" ]] || echo "WARNING ${meta.id}: contamination estimate unavailable - filtering without it" >&2
    fi
    G FilterMutectCalls -R ${fasta} -V merged.vcf.gz --ob-priors read-orientation-model.tar.gz \$EXTRA \\
        -O ${meta.id}.mutect2.filtered.vcf.gz ${params.filter_mutect_args}
    """
}

process MSISENSOR2 {
    tag "${meta.id}"
    label 'process_low'
    errorStrategy 'ignore'          // MSI is a bonus metric: never fail the run because of it

    input:
    tuple val(meta), path(bam), path(bai)
    path models

    output:
    tuple val(meta), path("${meta.id}.msi"), emit: msi

    script:
    """
    msisensor2 msi -M ${models} -t ${bam} -o ${meta.id}.msi -b ${task.cpus} -c ${params.msi_min_coverage}
    """
}
