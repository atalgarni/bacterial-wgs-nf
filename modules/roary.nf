// =============================================================================
//  Module: ROARY
//  Pangenome (core / accessory genes) from Prokka GFF3 files of all samples
//  that passed QC. Requires >= 2 genomes (enforced in main.nf).
// =============================================================================

process ROARY {
    tag "pangenome (${gffs.size()} genomes)"
    label 'process_high'

    publishDir "${params.outdir}/09_pangenome", mode: 'copy'

    input:
    path gffs

    output:
    path "roary", emit: results

    script:
    def aln = params.roary_core_alignment ? ('-e' + (params.roary_mafft ? ' --mafft' : '')) : ''
    """
    roary -f roary ${aln} -i ${params.roary_identity} -cd ${params.roary_core} -p ${task.cpus} *.gff
    """

    stub:
    """
    mkdir -p roary
    echo -e "Core genes\\t3000" > roary/summary_statistics.txt
    """
}
