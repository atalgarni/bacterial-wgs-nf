// =============================================================================
//  Modules: SNIPPY + SNIPPY_CORE
//  Variant calling of each assembly against the reference (contigs are
//  shredded into synthetic reads by snippy --ctgs), then a core-genome SNP
//  alignment across all samples. Equivalent to snippy-multi, but parallel
//  and portable (no absolute work-directory paths).
// =============================================================================

process SNIPPY {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/10_snippy", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)
    path reference

    output:
    path "${sample_id}", emit: dir

    script:
    """
    snippy --cpus ${task.cpus} --outdir ${sample_id} --prefix snps \\
        --ref ${reference} --ctgs ${contigs}
    """

    stub:
    """
    mkdir -p ${sample_id}; touch ${sample_id}/snps.tab
    """
}

process SNIPPY_CORE {
    tag "core (${dirs.size()} genomes)"
    label 'process_medium'

    publishDir "${params.outdir}/10_snippy/core", mode: 'copy'

    input:
    path dirs
    path reference

    output:
    path "core.*", emit: core

    script:
    """
    snippy-core --ref ${reference} --prefix core ${dirs}
    """

    stub:
    """
    touch core.aln core.full.aln core.tab core.txt core.vcf
    """
}
