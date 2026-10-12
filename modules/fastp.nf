// =============================================================================
//  Module: FASTP
//  Adapter detection/removal and quality trimming of paired-end reads.
// =============================================================================

process FASTP {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/02_fastp", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: reads
    path "${sample_id}_fastp.json", emit: json
    path "${sample_id}_fastp.html", emit: html

    script:
    """
    fastp \\
        -i ${r1} -I ${r2} \\
        -o ${sample_id}_R1.trimmed.fastq.gz -O ${sample_id}_R2.trimmed.fastq.gz \\
        --json ${sample_id}_fastp.json --html ${sample_id}_fastp.html \\
        --thread ${task.cpus} ${params.fastp_args}
    """

    stub:
    """
    cp ${r1} ${sample_id}_R1.trimmed.fastq.gz
    cp ${r2} ${sample_id}_R2.trimmed.fastq.gz
    echo '{}' > ${sample_id}_fastp.json; touch ${sample_id}_fastp.html
    """
}
