// =============================================================================
//  Module: SPECIES_CHECK
//  Purity / identity gate before assembly.
//  Bracken output is not ordered by abundance, so rows are sorted by
//  fraction_total_reads (column 7) and the most abundant species is taken.
//  PASS requires: top species fraction >= --min_species_fraction AND
//  (only if --target_species is set) the top species is the target species.
//  The top species is also emitted so that later steps (Prokka organism
//  label, Kaptive auto-selection) can use it per sample.
// =============================================================================

process SPECIES_CHECK {
    tag "${sample_id}"
    label 'process_single'

    input:
    tuple val(sample_id), path(bracken)

    output:
    tuple val(sample_id), env(STATUS), env(TOP), emit: status
    path "${sample_id}_species_check.tsv", emit: row

    script:
    def target = params.target_species ?: ''
    """
    awk -F'\\t' -v OFS='\\t' -v s="${sample_id}" -v target="${target}" -v min="${params.min_species_fraction}" '
        NR == 1 { next }
        {
            f = \$7 + 0
            if (f > best) { best = f; top = \$1 }
            if (target != "" && tolower(\$1) == tolower(target)) tf = f
        }
        END {
            ok_frac  = (best >= min + 0)
            ok_name  = (target == "" || tolower(top) == tolower(target))
            status   = (ok_frac && ok_name) ? "PASS" : "FAIL"
            reason   = status == "PASS" ? "-" : (!ok_name ? "top species is not the target" : "top species below threshold")
            print s, top, sprintf("%.4f", best), (target == "" ? "-" : target), (target == "" ? "-" : sprintf("%.4f", tf + 0)), status, reason > s "_species_check.tsv"
        }' ${bracken}
    STATUS=\$(cut -f6 ${sample_id}_species_check.tsv)
    TOP=\$(cut -f2 ${sample_id}_species_check.tsv)
    """

    stub:
    // Deterministic per-sample species so the test profile exercises both
    // Kaptive-supported and unsupported species.
    def opts = ['Acinetobacter baumannii', 'Klebsiella pneumoniae', 'Escherichia coli']
    def top  = opts[Math.abs(sample_id.hashCode()) % 3]
    """
    STATUS=PASS
    TOP="${top}"
    echo -e "${sample_id}\\t${top}\\t0.9900\\t-\\t-\\tPASS\\t-" > ${sample_id}_species_check.tsv
    """
}
