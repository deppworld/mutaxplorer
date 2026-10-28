// MutaXplorer - FASTQ pre-processing, alignment, duplicate / UMI handling, BQSR, BAM import

process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.trim_R1.fastq.gz"), path("${meta.id}.trim_R2.fastq.gz"), emit: reads
    path "${meta.id}.fastp.json", emit: json
    path "${meta.id}.fastp.html", emit: html

    script:
    // UMI at the start of a read is moved into the read name (":<UMI>"), where fgbio picks it up
    def umi = (params.umi && params.umi_loc) ? "--umi --umi_loc ${params.umi_loc} --umi_len ${params.umi_len} --umi_prefix '' ${params.umi_skip ? "--umi_skip ${params.umi_skip}" : ''}" : ''
    """
    fastp -i ${r1} -I ${r2} \\
        -o ${meta.id}.trim_R1.fastq.gz -O ${meta.id}.trim_R2.fastq.gz \\
        --detect_adapter_for_pe --length_required ${params.min_read_length} \\
        --thread ${task.cpus} ${umi} ${params.fastp_args} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
}

process ALIGN {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(r1), path(r2)
    tuple path(fasta), path(fai), path(dict)
    path index

    output:
    tuple val(meta), path("${meta.id}.raw.bam"), emit: bam     // query-name grouped

    script:
    def rg = "@RG\\tID:${meta.id}.${meta.lane ?: 1}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:ILLUMINA"
    """
    IDX=\$(ls ${index}/*.bwt.2bit.64 | head -1); IDX=\${IDX%.bwt.2bit.64}   # index prefix need not match the FASTA name
    bwa-mem2 mem -t ${task.cpus} -Y -R '${rg}' \$IDX ${r1} ${r2} \\
        | samtools view -@ 2 -b -o ${meta.id}.raw.bam -
    """
}

process MARK_DUPLICATES {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.md.bam"), path("${meta.id}.md.bam.bai"), emit: bam
    path "${meta.id}.markdup.txt", emit: metrics

    script:
    """
    samtools fixmate -@ ${task.cpus} -m ${bam} - \\
      | samtools sort -@ ${task.cpus} -m 1500M -T tmp_sort - \\
      | samtools markdup -@ ${task.cpus} -f ${meta.id}.markdup.txt -s - ${meta.id}.md.bam
    samtools index ${meta.id}.md.bam
    """
}

process UMI_CONSENSUS {
    // UMI read families -> one consensus read (fgbio), then re-alignment of the consensus reads
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(bam)
    tuple path(fasta), path(fai), path(dict)
    path index

    output:
    tuple val(meta), path("${meta.id}.umi.bam"), path("${meta.id}.umi.bam.bai"), emit: bam
    path "${meta.id}.umi_family_sizes.txt", emit: metrics

    script:
    def mem = (task.memory ? task.memory.toGiga() : 8) - 1
    def rg = "@RG\\tID:${meta.id}.umi\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:ILLUMINA"
    """
    samtools fixmate -m ${bam} fixmate.bam                 # mate tags (MQ/MC) needed by GroupReadsByUmi
    fgbio -Xmx${mem}g CopyUmiFromReadName -i fixmate.bam -o umi_tagged.bam
    fgbio -Xmx${mem}g GroupReadsByUmi -i umi_tagged.bam -o grouped.bam \\
        --strategy ${params.umi_strategy} --edits 1 --min-map-q 10 \\
        --family-size-histogram ${meta.id}.umi_family_sizes.txt
    fgbio -Xmx${mem}g CallMolecularConsensusReads -i grouped.bam -o consensus.unmapped.bam \\
        --min-reads ${params.umi_min_reads} --read-name-prefix ${meta.id} --threads ${task.cpus}
    samtools fastq -@ 2 -1 c_R1.fq.gz -2 c_R2.fq.gz -0 /dev/null -s /dev/null -n consensus.unmapped.bam
    IDX=\$(ls ${index}/*.bwt.2bit.64 | head -1); IDX=\${IDX%.bwt.2bit.64}
    bwa-mem2 mem -t ${task.cpus} -Y -R '${rg}' \$IDX c_R1.fq.gz c_R2.fq.gz \\
      | samtools sort -@ ${task.cpus} -m 1500M -T tmp_sort -o ${meta.id}.umi.bam -
    samtools index ${meta.id}.umi.bam
    rm -f fixmate.bam umi_tagged.bam grouped.bam consensus.unmapped.bam c_R*.fq.gz
    """
}

process BQSR {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)
    path known_sites        // list of VCFs (+ indexes)
    path intervals

    output:
    tuple val(meta), path("${meta.id}.recal.bam"), path("${meta.id}.recal.bam.bai"), emit: bam
    path "${meta.id}.recal.table", emit: table

    script:
    def known = known_sites.findAll { it.name.endsWith('.vcf.gz') || it.name.endsWith('.vcf') }.collect { "--known-sites ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" BaseRecalibrator \\
        -R ${fasta} -I ${bam} ${known} -L ${intervals} -O ${meta.id}.recal.table
    gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" ApplyBQSR \\
        -R ${fasta} -I ${bam} --bqsr-recal-file ${meta.id}.recal.table -O ${meta.id}.recal.bam
    mv ${meta.id}.recal.bai ${meta.id}.recal.bam.bai
    """
}

process PREPARE_BAM {
    // Imported BAM (e.g. DRAGEN / ICA): check contigs vs reference, guarantee RG/SM = sample id, index
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(bam)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.id}.input.bam"), path("${meta.id}.input.bam.bai"), emit: bam

    script:
    """
    samtools quickcheck ${bam} || { echo "ERROR ${meta.id}: BAM is truncated or corrupt" >&2; exit 1; }
    samtools view -H ${bam} | awk -F'\\t' '/^@SQ/{n="";l="";for(i=2;i<=NF;i++){if(\$i~/^SN:/)n=substr(\$i,4);if(\$i~/^LN:/)l=substr(\$i,4)}print n"\\t"l}' > bam_contigs.txt
    cut -f1,2 ${fai} > ref_contigs.txt
    if ! grep -x -F -f bam_contigs.txt ref_contigs.txt > shared.txt || [[ ! -s shared.txt ]]; then
        echo "ERROR ${meta.id}: BAM contigs do not match the reference (wrong genome build? hg19 vs hg38?)" >&2; exit 1
    fi
    if [[ "\$(grep -x -F -f ref_contigs.txt bam_contigs.txt)" != "\$(cat shared.txt)" ]]; then
        echo "ERROR ${meta.id}: contig ORDER differs from the reference - use the FASTA the BAM was aligned to" >&2; exit 1
    fi
    sm=\$(samtools view -H ${bam} | grep '^@RG' | tr '\\t' '\\n' | awk -F: '\$1=="SM"{print \$2}' | sort -u)
    if [[ "\$sm" == "${meta.id}" ]]; then
        ln -s ${bam} ${meta.id}.input.bam
    else
        # add/replace read groups so the sample name equals the samplesheet id (needed by Mutect2 -normal)
        samtools addreplacerg -@ ${task.cpus} -r "ID:${meta.id}" -r "SM:${meta.id}" -r "LB:${meta.id}" -r "PL:ILLUMINA" \\
            -m overwrite_all -o ${meta.id}.input.bam ${bam}
    fi
    samtools index -@ ${task.cpus} ${meta.id}.input.bam
    """
}
