// MutaXplorer - coverage and alignment QC

process MOSDEPTH {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)
    path regions
    val  thresholds

    output:
    tuple val(meta), path("${meta.id}.mosdepth.summary.txt"), path("${meta.id}.thresholds.bed.gz"), emit: cov
    path "${meta.id}.mosdepth.*", emit: mqc

    script:
    """
    mosdepth -t ${task.cpus} -n --mapq 20 --by ${regions} --thresholds ${thresholds} ${meta.id} ${bam}
    """
}

process SAMTOOLS_STATS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "${meta.id}.{stats,flagstat}", emit: mqc

    script:
    """
    samtools stats ${bam} > ${meta.id}.stats
    samtools flagstat ${bam} > ${meta.id}.flagstat
    """
}

process MULTIQC {
    label 'process_single'

    input:
    path files

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", optional: true

    script:
    """
    multiqc --force --title "MutaXplorer QC" --filename multiqc_report.html .
    """
}

process IDENTITY {
    // genotype every DNA sample at common SNPs and compare all pairs (sample swaps / mismatched pairs)
    label 'process_low'

    input:
    path bams
    path bais
    tuple path(fasta), path(fai), path(dict)
    tuple path(sites), path(sites_tbi)

    output:
    path "identity_genotypes.vcf.gz", emit: gtcheck

    script:
    """
    ls *.bam | sort > bams.txt
    bcftools mpileup -f ${fasta} -T ${sites} -a FORMAT/AD,FORMAT/DP -d 2000 -q 20 -Q 20 -b bams.txt -Ou 2>/dev/null \\
      | bcftools call -m -Oz -o identity_genotypes.vcf.gz
    """
}
