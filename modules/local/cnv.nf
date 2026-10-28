// MutaXplorer - copy-number analysis (CNVkit)

process CNVKIT_REFERENCE {
    // pooled reference from the normal samples in the run, or a flat reference when there are none
    tag "${normals ? normals.size() + ' normals' : 'flat'}"
    label 'process_medium'
    label 'cnvkit'

    input:
    path normals            // BAMs (may be empty)
    path normal_bais
    tuple path(fasta), path(fai), path(dict)
    path targets            // [] for WGS

    output:
    path "reference.cnn", emit: cnn

    script:
    def t = targets ? "--targets ${targets}" : ''
    def bams = normals instanceof List ? normals.join(' ') : normals.toString()
    """
    cnvkit.py batch --normal ${bams} ${t} --fasta ${fasta} --method ${params.cnv_method} \\
        --output-reference reference.cnn --output-dir ref_build -p ${task.cpus}
    """
}

process CNVKIT_CALL {
    tag "${meta.id}"
    label 'process_medium'
    label 'cnvkit'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("${meta.id}.call.cns"), path("${meta.id}.cnr"), emit: calls
    path "${meta.id}.cns", emit: cns
    path "${meta.id}-*.{png,pdf}", optional: true, emit: plots

    script:
    def purity = params.purity ? "--purity ${params.purity}" : ''
    """
    ln -s ${bam} ${meta.id}.bam && ln -s ${bai} ${meta.id}.bam.bai
    # use the R (with DNAcopy) that ships next to CNVkit, not whatever Rscript is first on PATH
    RS=\$(dirname \$(readlink -f \$(which cnvkit.py)))/Rscript; [[ -x \$RS ]] && RSARG="--rscript-path \$RS" || RSARG=""
    cnvkit.py batch ${meta.id}.bam -r ${reference} -d . -p ${task.cpus} --drop-low-coverage --scatter --diagram \$RSARG
    if [[ -n "${purity}" ]]; then
        cnvkit.py call ${meta.id}.cns ${purity} -o ${meta.id}.call.cns
    fi
    """
}
