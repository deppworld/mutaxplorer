// MutaXplorer - RNA fusion detection (STAR + Arriba)

process STAR_INDEX {
    tag "${fasta.name}"
    label 'process_high_memory'
    label 'arriba'

    input:
    tuple path(fasta), path(fai), path(dict)
    path gtf

    output:
    path "star", emit: index

    script:
    """
    GTF=${gtf}; [[ \$GTF == *.gz ]] && { zcat \$GTF > genes.gtf; GTF=genes.gtf; }
    LEN=\$(awk '{s+=\$2} END{print s}' ${fai})
    NB=\$(awk -v l=\$LEN 'BEGIN{n=int(log(l)/log(2)/2-1); print (n>14?14:n)}')
    mkdir star
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir star --genomeFastaFiles ${fasta} \\
         --sjdbGTFfile \$GTF --sjdbOverhang ${params.rna_read_length - 1} --genomeSAindexNbases \$NB
    """
}

process ARRIBA {
    tag "${meta.id}"
    label 'process_high'
    label 'arriba'

    input:
    tuple val(meta), path(r1), path(r2)
    tuple path(fasta), path(fai), path(dict)
    path gtf
    path star_index

    output:
    tuple val(meta), path("${meta.id}.arriba.tsv"), emit: fusions
    path "${meta.id}.arriba.discarded.tsv", emit: discarded

    script:
    def build = params.genome == 'hg19' ? 'hg19_hs37d5_GRCh37' : 'hg38_GRCh38'
    """
    DB=\$(dirname \$(which arriba))/../var/lib/arriba
    pick(){ [[ "\$1" == none ]] && return; [[ -n "\$1" ]] && { echo "\$1"; return; }; ls \$DB/\$2_${build}_*.\$3 2>/dev/null | tail -1; }
    BL=\$(pick "${params.arriba_blacklist ?: ''}" blacklist tsv.gz)
    KF=\$(pick "${params.arriba_known_fusions ?: ''}" known_fusions tsv.gz)
    PD=\$(pick "${params.arriba_protein_domains ?: ''}" protein_domains gff3)
    ARGS=""
    [[ -n "\$BL" ]] && ARGS="\$ARGS -b \$BL" || ARGS="\$ARGS -f blacklist"
    [[ -n "\$KF" ]] && ARGS="\$ARGS -k \$KF -t \$KF"
    [[ -n "\$PD" ]] && ARGS="\$ARGS -p \$PD"
    STAR --runThreadN ${task.cpus} --genomeDir ${star_index} --genomeLoad NoSharedMemory \\
        --readFilesIn ${r1} ${r2} --readFilesCommand zcat \\
        --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \\
        --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType WithinBAM HardClip \\
        --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 \\
        --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 \\
      | arriba -x /dev/stdin -o ${meta.id}.arriba.tsv -O ${meta.id}.arriba.discarded.tsv \\
          -a ${fasta} -g ${gtf} \$ARGS
    """
}
