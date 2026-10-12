// =============================================================================
//  Module: FASTANI
//  Average nucleotide identity of each assembly to the reference genome.
// =============================================================================

process FASTANI {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/07_fastani", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)
    path reference

    output:
    path "${sample_id}.ani.txt", emit: ani

    script:
    """
    fastANI -q ${contigs} -r ${reference} -o ${sample_id}.ani.txt --threads ${task.cpus}
    """

    stub:
    """
    echo -e "${sample_id}_contigs.fasta\\t${reference}\\t99.5\\t1200\\t1250" > ${sample_id}.ani.txt
    """
}
