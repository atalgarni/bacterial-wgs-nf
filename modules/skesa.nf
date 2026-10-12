// =============================================================================
//  Module: SKESA
//  De novo assembly of paired-end short reads (conservative k-mer extension).
// =============================================================================

process SKESA {
    tag "${sample_id}"
    label 'process_high'

    publishDir "${params.outdir}/06_assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta"), emit: contigs

    script:
    def mem_gb = task.memory ? task.memory.toGiga() : 16
    """
    skesa --reads ${r1},${r2} \\
        --cores ${task.cpus} --memory ${mem_gb} \\
        --min_contig ${params.min_contig} \\
        --contigs_out ${sample_id}_contigs.fasta
    """

    stub:
    """
    printf '>Contig_1\\nACGTACGTACGT\\n>Contig_2\\nGGCCAATT\\n' > ${sample_id}_contigs.fasta
    """
}
