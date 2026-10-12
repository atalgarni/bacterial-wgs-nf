// =============================================================================
//  Module: REMOVE_PHIX
//  Remove PhiX spike-in reads. Read pairs in which BOTH mates are unmapped to
//  PhiX are kept (samtools -f 12), so R1/R2 stay synchronised.
// =============================================================================

process REMOVE_PHIX {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/01_phix_removed", mode: 'copy', pattern: '*.{log,tsv}'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path phix_fasta

    output:
    tuple val(sample_id), path("${sample_id}_clean_R1.fastq.gz"), path("${sample_id}_clean_R2.fastq.gz"), emit: reads
    path "${sample_id}_phix.tsv",                                                                       emit: stats

    script:
    """
    bwa index ${phix_fasta} 2> /dev/null

    bwa mem -t ${task.cpus} ${phix_fasta} ${r1} ${r2} 2> ${sample_id}_phix.log \\
        | samtools view -b -f 12 -F 256 \\
        | samtools sort -n -@ ${task.cpus} -o ${sample_id}.unmapped.bam

    samtools fastq -@ ${task.cpus} \\
        -1 ${sample_id}_clean_R1.fastq.gz \\
        -2 ${sample_id}_clean_R2.fastq.gz \\
        -0 /dev/null -s /dev/null -n \\
        ${sample_id}.unmapped.bam

    RAW=\$(gunzip -c ${r1} | awk 'END{print NR/4}')
    KEPT=\$(gunzip -c ${sample_id}_clean_R1.fastq.gz | awk 'END{print NR/4}')
    echo -e "sample\\traw_pairs\\tkept_pairs\\tphix_pairs" > ${sample_id}_phix.tsv
    echo -e "${sample_id}\\t\${RAW}\\t\${KEPT}\\t\$((RAW-KEPT))" >> ${sample_id}_phix.tsv

    rm -f ${sample_id}.unmapped.bam
    """

    stub:
    """
    printf '@r/1\\nACGT\\n+\\nIIII\\n' | gzip > ${sample_id}_clean_R1.fastq.gz
    printf '@r/2\\nACGT\\n+\\nIIII\\n' | gzip > ${sample_id}_clean_R2.fastq.gz
    echo -e "sample\\traw_pairs\\tkept_pairs\\tphix_pairs\\n${sample_id}\\t1\\t1\\t0" > ${sample_id}_phix.tsv
    """
}
