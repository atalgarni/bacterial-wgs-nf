// =============================================================================
//  Module: FETCH_REFERENCE
//  Downloads a curated reference genome (assets/reference_genomes.tsv) from
//  NCBI with the Datasets CLI, VERIFIES that the assembly belongs to the
//  expected species, and keeps the chromosome(s) only (plasmids are dropped
//  so that core-SNP coordinates are chromosomal). The result is cached in
//  --reference_cache and reused by later runs (storeDir).
// =============================================================================

process FETCH_REFERENCE {
    tag "${key}"
    label 'process_single'

    storeDir "${params.reference_cache}/${key}"

    input:
    tuple val(key), val(species), val(strain), val(accession)

    output:
    tuple val(key), val(species), path("${key}.fasta"), emit: fasta
    path "${key}.info.tsv",                               emit: info

    script:
    """
    datasets download genome accession ${accession} \\
        --include genome,seq-report --filename ref.zip
    unzip -q -o ref.zip -d ref

    DATA=\$(ls -d ref/ncbi_dataset/data/GC*_* | head -1)
    FNA=\$(ls \$DATA/*.fna | head -1)

    # 1. Species check: the downloaded assembly must be the expected species
    ORG=\$(grep -o '"organismName":"[^"]*"' ref/ncbi_dataset/data/assembly_data_report.jsonl | head -1 | cut -d'"' -f4)
    if ! echo "\$ORG" | grep -qi "^${species}"; then
        echo "ERROR: ${accession} is '\$ORG', expected '${species}'. Fix assets/reference_genomes.tsv." >&2
        exit 1
    fi

    # 2. Keep chromosome(s); fall back to all sequences if roles are absent
    grep '"role":"assembled-molecule"' \$DATA/sequence_report.jsonl 2>/dev/null \\
        | grep -i '"assignedMoleculeLocationType":"Chromosome"' \\
        | grep -o '"refseqAccession":"[^"]*"' | cut -d'"' -f4 > chrom.ids || true
    if [ -s chrom.ids ]; then
        awk 'NR==FNR { keep[\$1] = 1; next }
             /^>/    { id = substr(\$1, 2); p = (id in keep) }
             p' chrom.ids \$FNA > ${key}.fasta
        CONTENT="chromosome(s): \$(paste -sd, chrom.ids)"
    else
        cp \$FNA ${key}.fasta
        CONTENT="all sequences (no chromosome roles in report)"
    fi

    LEN=\$(grep -v '^>' ${key}.fasta | tr -d '\\n' | wc -c)
    [ "\$LEN" -gt 100000 ] || { echo "ERROR: reference ${key} is only \$LEN bp" >&2; exit 1; }

    printf 'key\\tspecies\\tstrain\\taccession\\torganism_in_ncbi\\tcontent\\tlength_bp\\n' > ${key}.info.tsv
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "${key}" "${species}" "${strain}" "${accession}" "\$ORG" "\$CONTENT" "\$LEN" >> ${key}.info.tsv
    rm -rf ref ref.zip
    """

    stub:
    """
    printf '>STUB_${key}\\nACGTACGTACGTACGTACGT\\n' > ${key}.fasta
    printf 'key\\tspecies\\tstrain\\taccession\\torganism_in_ncbi\\tcontent\\tlength_bp\\n${key}\\t${species}\\t${strain}\\t${accession}\\t${species}\\tstub\\t20\\n' > ${key}.info.tsv
    """
}
