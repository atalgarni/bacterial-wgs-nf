// =============================================================================
//  Module: KRAKEN2
//  k-mer based taxonomic classification of trimmed reads (report only).
// =============================================================================

process KRAKEN2 {
    tag "${sample_id}"
    label 'process_high'

    publishDir "${params.outdir}/03_kraken2", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path db

    output:
    tuple val(sample_id), path("${sample_id}.kraken2.report"), emit: report

    script:
    """
    kraken2 --db ${db} --threads ${task.cpus} \\
        --report ${sample_id}.kraken2.report \\
        --paired ${r1} ${r2} > /dev/null
    """

    stub:
    """
    touch ${sample_id}.kraken2.report
    """
}
