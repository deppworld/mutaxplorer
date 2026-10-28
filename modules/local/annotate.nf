// MutaXplorer - normalisation and annotation

process NORMALIZE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.id}.norm.vcf.gz"), path("${meta.id}.norm.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools norm -m -any -f ${fasta} -Oz -o ${meta.id}.norm.vcf.gz ${vcf}
    bcftools index -t ${meta.id}.norm.vcf.gz
    """
}

process VEP {
    tag "${meta.id}"
    label 'process_medium'
    label 'vep'

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)
    path cache
    val  assembly

    output:
    tuple val(meta), path("${meta.id}.annotated.vcf.gz"), emit: vcf

    script:
    """
    bcftools view -f PASS -Oz -o pass.vcf.gz ${vcf} 2>/dev/null || zcat ${vcf} | awk '/^#/||\$7=="PASS"' | bgzip > pass.vcf.gz
    vep --offline --cache --dir_cache ${cache} --assembly ${assembly} --fasta ${fasta} \\
        --input_file pass.vcf.gz --output_file ${meta.id}.annotated.vcf.gz --vcf --compress_output bgzip \\
        --everything --pick --pick_order mane_select,canonical,appris,tsl,biotype,ccds,rank,length \\
        --fork ${task.cpus} --force_overwrite --no_stats ${params.vep_args}
    """
}
