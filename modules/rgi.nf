// =============================================================================
//  Modules: RGI_LOAD + RGI_MAIN
//  Antimicrobial resistance genes against CARD (Resistance Gene Identifier).
//  The CARD JSON is loaded once into a local database that every sample uses.
// =============================================================================

process RGI_LOAD {
    tag "CARD"
    label 'process_single'

    input:
    path card_json

    output:
    path "localDB", emit: db

    script:
    """
    rgi load --card_json ${card_json} --local
    """

    stub:
    """
    mkdir -p localDB
    """
}

process RGI_MAIN {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/11_rgi", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)
    path localDB

    output:
    path "${sample_id}.txt",  emit: txt
    path "${sample_id}.json", emit: json

    script:
    """
    rgi main --input_sequence ${contigs} --output_file ${sample_id} \\
        --input_type contig --alignment_tool DIAMOND --data wgs \\
        --local --clean --num_threads ${task.cpus}
    """

    stub:
    """
    touch ${sample_id}.txt; echo '{}' > ${sample_id}.json
    """
}
