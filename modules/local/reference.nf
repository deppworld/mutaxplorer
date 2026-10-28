// MutaXplorer - reference preparation and interval handling

process PREPARE_REFERENCE {
    tag "${fasta.name}"
    label 'process_low'

    input:
    path fasta

    output:
    tuple path(fasta), path("${fasta}.fai"), path("${fasta.baseName}.dict"), emit: ref

    script:
    """
    [[ -s ${fasta}.fai ]] || samtools faidx ${fasta}
    [[ -s ${fasta.baseName}.dict ]] || samtools dict ${fasta} -o ${fasta.baseName}.dict
    """
}

process BWA_INDEX {
    tag "${fasta.name}"
    label 'process_high_memory'

    input:
    path fasta

    output:
    path "bwamem2", emit: index

    script:
    """
    mkdir bwamem2
    ln -s ../${fasta} bwamem2/${fasta}
    bwa-mem2 index bwamem2/${fasta}
    """
}

process MAKE_INTERVALS {
    tag "${params.assay}"
    label 'process_single'

    input:
    path fai
    path targets      // [] for WGS
    path exclude      // [] or BED
    val  chunks

    output:
    path "intervals.all.bed",       emit: all
    path "intervals.chunk_*.bed",   emit: chunks

    script:
    def bed = targets ? "--bed ${targets} --pad ${params.interval_padding}" : ''
    def ex  = exclude ? "--exclude ${exclude}" : ''
    """
    split_intervals.py --fai ${fai} ${bed} ${ex} --chunks ${chunks} --prefix intervals
    """
}

process GENE_BED {
    tag "${gtf.name}"
    label 'process_single'

    input:
    path gtf
    path fai

    output:
    path "genes.bed", emit: bed

    script:
    """
    gtf_to_genes.py --gtf ${gtf} --fai ${fai} --out genes.bed
    """
}
