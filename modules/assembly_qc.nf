// =============================================================================
//  Module: ASSEMBLY_QC
//  Assembly metrics (contigs, total length, largest, N50, GC) and the
//  contig-count gate (--max_contigs). POSIX awk only.
// =============================================================================

process ASSEMBLY_QC {
    tag "${sample_id}"
    label 'process_single'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), env(STATUS), emit: status
    path "${sample_id}_assembly_qc.tsv", emit: row

    script:
    """
    awk '/^>/ { if (s != "") print length(s) "\\t" gsub(/[GCgc]/, "", s); s = ""; next }
         { s = s \$0 }
         END { if (s != "") print length(s) "\\t" gsub(/[GCgc]/, "", s) }' ${contigs} \\
      | sort -k1,1nr \\
      | awk -F'\\t' -v OFS='\\t' -v s="${sample_id}" -v max="${params.max_contigs}" '
            { l[NR] = \$1; t += \$1; gc += \$2 }
            END {
                r = 0; n50 = 0
                for (i = 1; i <= NR; i++) { r += l[i]; if (!n50 && r >= t / 2) n50 = l[i] }
                status = (NR > 0 && NR <= max + 0) ? "PASS" : "FAIL"
                print s, NR, t, (NR ? l[1] : 0), n50, (t ? sprintf("%.2f", gc / t * 100) : "0"), status
            }' > ${sample_id}_assembly_qc.tsv
    STATUS=\$(cut -f7 ${sample_id}_assembly_qc.tsv)
    """

    stub:
    """
    STATUS=PASS
    echo -e "${sample_id}\\t2\\t20\\t12\\t12\\t50.00\\tPASS" > ${sample_id}_assembly_qc.tsv
    """
}
