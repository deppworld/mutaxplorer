// MutaXplorer - germline SNV/indel calling (GATK4 HaplotypeCaller + hard filters)

process HAPLOTYPECALLER {
    tag "${meta.id}:${interval.baseName.tokenize('.').last()}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), path(interval)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: calls

    script:
    prefix = "${meta.id}.${interval.baseName.tokenize('.').last()}"
    """
    gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" HaplotypeCaller \\
        -R ${fasta} -I ${bam} -L ${interval} \\
        --native-pair-hmm-threads ${task.cpus} \\
        -O ${prefix}.vcf.gz ${params.haplotypecaller_args}
    """
}

process GERMLINE_FILTER {
    // merge interval chunks, then GATK-recommended hard filters (per variant type)
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(vcfs), path(tbis)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.id}.germline.filtered.vcf.gz"), path("${meta.id}.germline.filtered.vcf.gz.tbi"), emit: vcf

    script:
    def ins = (vcfs instanceof List ? vcfs : [vcfs]).collect { "-I ${it}" }.join(' ')
    """
    G() { gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" "\$@"; }
    G MergeVcfs ${ins} -O merged.vcf.gz
    G VariantFiltration -R ${fasta} -V merged.vcf.gz -O ${meta.id}.germline.filtered.vcf.gz \\
        --filter-name SNP_QD    --filter-expression "vc.isSNP() && QD < 2.0" \\
        --filter-name SNP_FS    --filter-expression "vc.isSNP() && FS > 60.0" \\
        --filter-name SNP_MQ    --filter-expression "vc.isSNP() && MQ < 40.0" \\
        --filter-name SNP_SOR   --filter-expression "vc.isSNP() && SOR > 3.0" \\
        --filter-name SNP_MQRS  --filter-expression "vc.isSNP() && MQRankSum < -12.5" \\
        --filter-name SNP_RPRS  --filter-expression "vc.isSNP() && ReadPosRankSum < -8.0" \\
        --filter-name INDEL_QD  --filter-expression "!vc.isSNP() && QD < 2.0" \\
        --filter-name INDEL_FS  --filter-expression "!vc.isSNP() && FS > 200.0" \\
        --filter-name INDEL_SOR --filter-expression "!vc.isSNP() && SOR > 10.0" \\
        --filter-name INDEL_RPRS --filter-expression "!vc.isSNP() && ReadPosRankSum < -20.0" \\
        --genotype-filter-expression "GQ < 20" --genotype-filter-name LowGQ
    """
}
