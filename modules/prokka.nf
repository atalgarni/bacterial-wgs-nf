// =============================================================================
//  Module: PROKKA
//  Gene prediction and functional annotation of each assembly.
//  Organism label (--genus/--species) is taken, in order of priority, from
//  --prokka_genus/--prokka_species, else the sample's species (Bracken top
//  species or --species). Unknown species → Prokka's generic defaults.
// =============================================================================

process PROKKA {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/08_prokka", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs), val(organism)
    path proteins

    output:
    tuple val(sample_id), path("${sample_id}"), emit: dir
    path "${sample_id}/${sample_id}.gff",        emit: gff

    script:
    // "Acinetobacter baumannii" → genus/species; brackets and odd characters dropped
    def words   = (organism ?: '').replaceAll(/[^A-Za-z0-9_. -]/, '').tokenize(' ')
    def g       = params.prokka_genus   ?: (words.size() >= 1 && organism != 'unknown' ? words[0] : null)
    def sp      = params.prokka_species ?: (words.size() >= 2 && organism != 'unknown' ? words[1] : null)
    def genus   = g  ? "--genus ${g}" + (params.prokka_usegenus ? ' --usegenus' : '') : ''
    def species = sp ? "--species ${sp}" : ''
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
