// MutaXplorer - cohort aggregation: SNV/indel, TMB, MSI, CNV, SV, fusions, signatures, identity -> report

process COHORT_REPORT {
    tag "${mode}"
    label 'process_low'

    input:
    val  mode               // somatic | germline
    path files              // all per-sample result files (flat)
    path samples            // samples.tsv: id, patient, status, type
    tuple path(fasta), path(fai), path(dict)
    path genes              // genes.bed or []
    path intervals          // calling intervals
    path known_fusions      // [] or file
    path signatures         // [] or COSMIC-style SBS matrix

    output:
    path "${mode}/cohort/*", emit: tables
    path "${mode}/report/*", emit: report

    script:
    def vaf = mode == 'somatic' ? params.min_vaf : params.germline_min_vaf
    def tmb = params.tmb_region_mb ? "--tmb-region-mb ${params.tmb_region_mb}" : ''
    def tg  = params.assay == 'wgs' ? '' : "--targets ${intervals}"
    def sub = "${params.assay.toUpperCase()} | ${params.genome} | ${mode == 'somatic' ? 'Mutect2, CNVkit, Manta, Arriba' : 'HaplotypeCaller, CNVkit, Manta'} + VEP"
    """
    mkdir -p calls ${mode}/cohort ${mode}/report
    for f in ${files}; do ln -sf \$(readlink -f \$f) calls/\$(basename \$f); done
    summarize_structural.py --calls-dir calls --tables ${mode}/cohort --report ${mode}/report --mode ${mode} \\
        --samples ${samples} ${genes ? "--gene-bed ${genes}" : ''} ${tg} --fasta ${fasta} \\
        ${known_fusions ? "--known-fusions ${known_fusions}" : ''} ${signatures ? "--signatures ${signatures}" : ''} \\
        --amp-cn ${params.cnv_amp_cn} --sv-min-support ${params.sv_min_support} \\
        --signature-min-snvs ${params.signature_min_snvs} --min-dp ${params.min_dp}
    aggregate_cohort.py --mode ${mode} --calls-dir calls --out ${mode}/cohort \\
        --min-dp ${params.min_dp} --min-alt ${params.min_alt_reads} --min-vaf ${vaf} \\
        --max-pop-af ${params.max_pop_af} --msi-high ${params.msi_high_cutoff} \\
        --msi-min-sites ${params.msi_min_sites} ${tmb}
    visualize_cohort.py --mode ${mode} --cohort ${mode}/cohort --out ${mode}/report \\
        --subtitle "${sub}" --top-genes ${params.top_n_genes} --top-variants ${params.top_n_variants} \\
        --msi-high ${params.msi_high_cutoff} --min-dp ${params.min_dp} --min-alt ${params.min_alt_reads} \\
        --min-vaf ${vaf} --max-pop-af ${params.max_pop_af} \\
        --events ${mode}/cohort/oncoplot_events.tsv --extra-html ${mode}/report/sections.html
    rm -f ${mode}/report/sections.html
    """
}
