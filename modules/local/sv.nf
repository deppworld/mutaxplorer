// MutaXplorer - structural variants (Manta): somatic (tumor-only / tumor-normal) and germline

process MANTA {
    tag "${meta.id}"
    label 'process_medium'
    label 'manta'

    input:
    tuple val(meta), path(tbam), path(tbai), path(nbam), path(nbai)
    tuple path(fasta), path(fai), path(dict)
    path regions
    val  mode               // somatic | germline

    output:
    tuple val(meta), path("${meta.id}.manta.vcf.gz"), path("${meta.id}.manta.vcf.gz.tbi"), emit: vcf

    script:
    def exome = params.assay == 'wgs' ? '' : '--exome'
    def input = mode == 'germline' ? "--bam ${tbam}" : "--tumorBam ${tbam}" + (nbam ? " --normalBam ${nbam}" : '')
    def result = mode == 'germline' ? 'diploidSV' : (nbam ? 'somaticSV' : 'tumorSV')
    """
    bgzip -c ${regions} > call_regions.bed.gz && tabix -p bed call_regions.bed.gz
    configManta.py ${input} --referenceFasta ${fasta} --runDir manta ${exome} --callRegions call_regions.bed.gz
    manta/runWorkflow.py -j ${task.cpus} -g ${task.memory.toGiga()}
    cp manta/results/variants/${result}.vcf.gz ${meta.id}.manta.vcf.gz
    cp manta/results/variants/${result}.vcf.gz.tbi ${meta.id}.manta.vcf.gz.tbi
    """
}
