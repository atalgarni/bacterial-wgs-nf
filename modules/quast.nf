// =============================================================================
//  Module: QUAST
//  Assembly quality report across all passing assemblies (reference-based
//  metrics when --reference is given).
// =============================================================================

process QUAST {
    tag "quast"
    label 'process_medium'

    publishDir "${params.outdir}/12_quast", mode: 'copy'

    input:
    path contigs
    path reference

    output:
    path "quast", emit: results

    script:
    def ref = reference.name != 'NO_FILE' ? "-r ${reference}" : ''
    """
    # QUAST runs jobs through joblib's multiprocessing backend and relies on module
    # globals reaching the workers. Python 3.8+ on macOS starts workers with 'spawn',
    # which loses them (ablab/quast#175); use 'fork' (the Linux default) instead.
    export JOBLIB_START_METHOD=fork
    quast.py ${contigs} ${ref} -o quast --threads ${task.cpus}
    """

    stub:
    """
    mkdir -p quast; touch quast/report.tsv
    """
}
