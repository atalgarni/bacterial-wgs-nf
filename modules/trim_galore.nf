// =============================================================================
//  Module: TRIM_GALORE
//  Cutadapt-based adapter and quality trimming of paired-end reads, with
//  FastQC on the trimmed reads (default trimmer; see --trimmer).
// =============================================================================

process TRIM_GALORE {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/02_trim_galore", mode: 'copy', pattern: '*.{txt,html,zip}'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_val_1.fq.gz"), path("${sample_id}_val_2.fq.gz"), emit: reads
    path "*_trimming_report.txt",                                                          emit: report
    path "*_fastqc.{html,zip}",                                                            emit: fastqc, optional: true

    script:
    // Trim Galore's --cores N spawns roughly 3N+3 threads (pigz + cutadapt);
    // keep N in 1..8 and leave headroom, as recommended in the Trim Galore docs.
    def cores = Math.min(8, Math.max(1, (task.cpus as int) - 4))
    def fastqc = params.trim_galore_fastqc ? "--fastqc" : ""
    """
    trim_galore \\
        --paired \\
        --cores ${cores} \\
        --quality ${params.trim_quality} \\
        --length ${params.trim_min_length} \\
        --basename ${sample_id} \\
        --gzip \\
        ${fastqc} \\
        ${params.trim_galore_args} \\
        ${r1} ${r2}
    """

    stub:
    """
    cp ${r1} ${sample_id}_val_1.fq.gz
    cp ${r2} ${sample_id}_val_2.fq.gz
    touch ${r1}_trimming_report.txt ${r2}_trimming_report.txt
    touch ${sample_id}_val_1_fastqc.html ${sample_id}_val_1_fastqc.zip
    """
}
