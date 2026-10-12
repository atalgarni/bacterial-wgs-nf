// =============================================================================
//  Modules: QC_SUMMARY + MULTIQC
// =============================================================================

process QC_SUMMARY {
    label 'process_single'

    publishDir "${params.outdir}/05_qc_gates", mode: 'copy'

    input:
    path species_rows
    path assembly_rows

    output:
    path "species_check.tsv"
    path "assembly_qc.tsv"

    script:
    """
    echo -e "sample\\ttop_species\\ttop_fraction\\ttarget_species\\ttarget_fraction\\tstatus\\treason" > species_check.tsv
    cat ${species_rows} | sort >> species_check.tsv
    echo -e "sample\\tcontigs\\ttotal_bp\\tlargest_contig\\tn50\\tgc_pct\\tstatus" > assembly_qc.tsv
    cat ${assembly_rows} | sort >> assembly_qc.tsv
    """
}

process MULTIQC {
    label 'process_single'

    publishDir "${params.outdir}/13_multiqc", mode: 'copy'

    input:
    path reports

    output:
    path "multiqc_report.html", optional: true
    path "multiqc_data",        optional: true

    script:
    """
    multiqc . --force
    """

    stub:
    """
    touch multiqc_report.html
    """
}
