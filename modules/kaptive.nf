// =============================================================================
//  Module: KAPTIVE
//  Surface-polysaccharide locus typing with Kaptive 3. SPECIES-RESTRICTED:
//  Kaptive databases exist only for
//    · Klebsiella pneumoniae species complex (KpSC): kpsc_k (K locus),
//                                                   kpsc_o (O locus)
//    · Acinetobacter baumannii:                     ab_k  (K locus, KL),
//                                                   ab_o  (outer-core locus, OCL)
//  One task per (sample, database). Which databases a sample receives is
//  decided in main.nf (--kaptive auto | off | comma-separated keywords).
// =============================================================================

process KAPTIVE {
    tag "${sample_id}:${db}"
    label 'process_low'

    publishDir "${params.outdir}/14_kaptive/per_sample", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs), val(db)

    output:
    tuple val(db), path("${sample_id}_${db}.tsv"), emit: tsv

    script:
    """
    kaptive assembly ${db} ${contigs} \\
        -o ${sample_id}_${db}.tsv \\
        -t ${task.cpus}

    # Kaptive writes only the header when no locus genes align (wrong species
    # or locus missing from the assembly). Record the sample as untypeable so
    # it is not silently absent from the summary table.
    if [ \$(wc -l < ${sample_id}_${db}.tsv) -lt 2 ]; then
        awk -F'\\t' -v OFS='\\t' -v s="${sample_id}" 'NR == 1 {
            \$0 = ""; \$1 = s; \$2 = "none"; \$3 = "none"; \$4 = "Untypeable"
            \$5 = "no locus gene alignments (check species)"; NF = 19; print
        }' ${sample_id}_${db}.tsv > untyped.row
        cat untyped.row >> ${sample_id}_${db}.tsv && rm untyped.row
    fi
    """

    stub:
    def locus = [ab_k: 'KL2', ab_o: 'OCL1', kpsc_k: 'KL64', kpsc_o: 'O1/O2v1'][db] ?: 'unknown'
    """
    printf 'Assembly\\tBest match locus\\tBest match type\\tMatch confidence\\tProblems\\tIdentity\\n' >  ${sample_id}_${db}.tsv
    printf '${sample_id}_contigs\\t${locus}\\tunknown\\tTypeable\\t\\t99.50%%\\n'                    >> ${sample_id}_${db}.tsv
    """
}
