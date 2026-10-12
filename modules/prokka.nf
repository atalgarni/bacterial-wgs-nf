// =============================================================================
//  Module: PROKKA
//  Gene prediction and functional annotation of each assembly.
// =============================================================================

process PROKKA {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/08_prokka", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)
    path proteins

    output:
    tuple val(sample_id), path("${sample_id}"), emit: dir
    path "${sample_id}/${sample_id}.gff",        emit: gff

    script:
    def genus      = params.prokka_genus   ? "--genus ${params.prokka_genus} --usegenus" : ''
    def species    = params.prokka_species ? "--species ${params.prokka_species}"        : ''
    def prot       = proteins.name != 'NO_FILE' ? "--proteins ${proteins}"              : ''
    """
    prokka --outdir ${sample_id} --force --prefix ${sample_id} \\
        --addgenes --centre ${params.prokka_centre} --compliant \\
        ${genus} ${species} --strain ${sample_id} \\
        ${prot} --evalue ${params.prokka_evalue} --rfam \\
        --cpus ${task.cpus} \\
        ${contigs}
    """

    stub:
    """
    mkdir -p ${sample_id}
    printf '##gff-version 3\\n' > ${sample_id}/${sample_id}.gff
    """
}
