// =============================================================================
//  Module: BRACKEN
//  Bayesian re-estimation of species abundance from the Kraken2 report.
// =============================================================================

process BRACKEN {
    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/04_bracken", mode: 'copy'

    input:
    tuple val(sample_id), path(kraken_report)
    path db

    output:
    tuple val(sample_id), path("${sample_id}.bracken"), emit: abundance
    path "${sample_id}_bracken.report",                 emit: report

    script:
    """
    bracken -d ${db} -i ${kraken_report} \\
        -o ${sample_id}.bracken -w ${sample_id}_bracken.report \\
        -r ${params.read_len} -l ${params.bracken_level}
    """

    stub:
    """
    printf 'name\\ttaxonomy_id\\ttaxonomy_lvl\\tkraken_assigned_reads\\tadded_reads\\tnew_est_reads\\tfraction_total_reads\\n' > ${sample_id}.bracken
    printf 'Escherichia coli\\t562\\tS\\t10\\t0\\t10\\t0.01000\\n' >> ${sample_id}.bracken
    printf '${params.target_species ?: "Acinetobacter baumannii"}\\t470\\tS\\t990\\t0\\t990\\t0.99000\\n' >> ${sample_id}.bracken
    touch ${sample_id}_bracken.report
    """
}
