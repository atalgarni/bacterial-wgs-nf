#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
//  bacterial-wgs-nf
//  Illumina paired-end bacterial isolate WGS:
//    PhiX removal → Trim Galore (or fastp) → Kraken2/Bracken → species gate → SKESA →
//    assembly QC gate → FastANI · Prokka · RGI · Snippy → Roary · snippy-core
//    · Kaptive (A. baumannii / K. pneumoniae SC only) · QUAST · MultiQC
//  Species-agnostic: the organism is taken per sample from Bracken (or
//  --species); only Kaptive is restricted to the species it has databases for.
// ─────────────────────────────────────────────────────────────────────────────

include { REMOVE_PHIX           } from './modules/remove_phix'
include { TRIM_GALORE           } from './modules/trim_galore'
include { FASTP                 } from './modules/fastp'
include { KRAKEN2               } from './modules/kraken2'
include { BRACKEN               } from './modules/bracken'
include { SPECIES_CHECK         } from './modules/species_check'
include { SKESA                 } from './modules/skesa'
include { ASSEMBLY_QC           } from './modules/assembly_qc'
include { FASTANI               } from './modules/fastani'
include { PROKKA                } from './modules/prokka'
include { ROARY                 } from './modules/roary'
include { SNIPPY; SNIPPY_CORE   } from './modules/snippy'
include { RGI_LOAD; RGI_MAIN    } from './modules/rgi'
include { QUAST                 } from './modules/quast'
include { KAPTIVE               } from './modules/kaptive'
include { QC_SUMMARY; MULTIQC   } from './modules/summary'

// ── Kaptive: the only species-restricted step ────────────────────────────────

KAPTIVE_DBS = ['ab_k', 'ab_o', 'kpsc_k', 'kpsc_o']

KPSC = ['klebsiella pneumoniae', 'klebsiella quasipneumoniae', 'klebsiella variicola',
        'klebsiella quasivariicola', 'klebsiella africana']

def kaptive_dbs_for(String organism) {
    def s = (organism ?: '').toLowerCase().trim()
    if (s == 'acinetobacter baumannii') return ['ab_k', 'ab_o']
    if (KPSC.any { s == it || s.startsWith(it + ' ') }) return ['kpsc_k', 'kpsc_o']
    return []
}

// ── Parameter validation ──────────────────────────────────────────────────────

def check_path(name, path, type) {
    if (!path) return
    def f = file(path)
    if (!f.exists()) error "ERROR: --${name} not found: ${path}"
    if (type == 'dir' && !f.isDirectory()) error "ERROR: --${name} must be a directory: ${path}"
}

def validate_params() {
    if (!params.reads && !params.assemblies)
        error "ERROR: provide --reads \"dir/*_R{1,2}*.fastq.gz\" or --assemblies \"dir/*.fasta\""
    if (params.reads && !params.kraken_db && !params.skip_species_check)
        error "ERROR: --kraken_db is required (or use --skip_species_check)"
    if (!(params.trimmer in ['trim_galore', 'fastp']))
        error "ERROR: --trimmer must be 'trim_galore' or 'fastp' (got '${params.trimmer}')"
    if (!(params.kaptive in ['auto', 'off'])) {
        def bad = params.kaptive.toString().tokenize(',')*.trim().findAll { !(it in KAPTIVE_DBS) }
        if (bad) error "ERROR: --kaptive must be 'auto', 'off' or a comma-separated list of ${KAPTIVE_DBS.join(', ')} (got: ${bad.join(', ')})"
    }
    if (workflow.stubRun) return
    check_path('kraken_db',  params.kraken_db,  'dir')
    check_path('phix',       params.phix,       'file')
    check_path('reference',  params.reference,  'file')
    check_path('card_json',  params.card_json,  'file')
    check_path('prokka_proteins', params.prokka_proteins, 'file')
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow {
    validate_params()

    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║  bacterial-wgs-nf  v${workflow.manifest.version}                                  ║
    ╚══════════════════════════════════════════════════════════╝
    Input           : ${params.reads ?: params.assemblies + ' (assemblies)'}
    Output dir      : ${params.outdir}
    PhiX removal    : ${params.reads && params.phix ? params.phix : 'off'}
    Read trimming   : ${params.trimmer}
    Species gate    : ${params.skip_species_check ? 'off' : "${params.target_species ?: 'any species (purity only)'} ≥ ${params.min_species_fraction}"}
    Species label   : ${params.reads && !params.skip_species_check ? 'per sample from Bracken' : (params.species ?: 'unknown (set --species)')}
    Contig gate     : ≤ ${params.max_contigs} contigs
    Reference       : ${params.reference ?: 'none (FastANI and Snippy skipped)'}
    CARD JSON       : ${params.card_json ?: 'none (RGI skipped)'}
    Pangenome       : ${params.skip_pangenome ? 'off' : 'Roary'}
    Kaptive         : ${params.kaptive == 'auto' ? 'auto (A. baumannii / K. pneumoniae SC samples only)' : params.kaptive}
    Profile         : ${workflow.profile}
    ══════════════════════════════════════════════════════════
    """.stripIndent()

    no_file    = file("${projectDir}/assets/NO_FILE")
    ch_ref     = params.reference ? Channel.value(file(params.reference)) : null
    ch_qc_rows = Channel.empty()
    ch_sp_rows = Channel.empty()
    ch_mqc     = Channel.empty()
    ch_species = null

    if (params.reads) {

        ch_reads = Channel
            .fromFilePairs(params.reads, flat: true, checkIfExists: true)
            .ifEmpty { error "ERROR: no read pairs match --reads ${params.reads}" }

        // 1. PhiX removal (optional)
        if (params.phix) {
            REMOVE_PHIX(ch_reads, Channel.value(file(params.phix)))
            ch_reads = REMOVE_PHIX.out.reads
        }

        // 2. Trimming
        if (params.trimmer == 'fastp') {
            FASTP(ch_reads)
            ch_trimmed = FASTP.out.reads
            ch_mqc = ch_mqc.mix(FASTP.out.json)
        } else {
            TRIM_GALORE(ch_reads)
            ch_trimmed = TRIM_GALORE.out.reads
            ch_mqc = ch_mqc.mix(TRIM_GALORE.out.report.flatten(), TRIM_GALORE.out.fastqc.flatten())
        }

        // 3–4. Taxonomy + species gate
        if (params.skip_species_check) {
            ch_for_assembly = ch_trimmed
        } else {
            ch_db = Channel.value(file(params.kraken_db))
            KRAKEN2(ch_trimmed, ch_db)
            BRACKEN(KRAKEN2.out.report, ch_db)
            SPECIES_CHECK(BRACKEN.out.abundance)
            ch_sp_rows = SPECIES_CHECK.out.row
            ch_species = SPECIES_CHECK.out.status.map { id, status, top -> tuple(id, top ?: 'unknown') }
            ch_mqc = ch_mqc.mix(KRAKEN2.out.report.map { it[1] }, BRACKEN.out.report)

            SPECIES_CHECK.out.status
                .filter { id, status, top -> status != 'PASS' }
                .subscribe { id, status, top -> log.warn "Species gate: ${id} FAILED — not assembled (see 05_qc_gates/species_check.tsv)" }

            ch_for_assembly = ch_trimmed
                .join(SPECIES_CHECK.out.status.map { id, status, top -> tuple(id, status) })
                .filter { id, r1, r2, status -> status == 'PASS' }
                .map    { id, r1, r2, status -> tuple(id, r1, r2) }
        }

        // 5. Assembly
        SKESA(ch_for_assembly)
        ch_contigs = SKESA.out.contigs

    } else {
        // Start from existing assemblies
        ch_contigs = Channel
            .fromPath(params.assemblies, checkIfExists: true)
            .map { f -> tuple(f.name.replaceAll(/(_contigs)?\.(fasta|fa|fna)$/, ''), f) }
    }

    // 6. Assembly QC gate
    ASSEMBLY_QC(ch_contigs)
    ch_qc_rows = ASSEMBLY_QC.out.row

    ASSEMBLY_QC.out.status
        .filter { id, status -> status != 'PASS' }
        .subscribe { id, status -> log.warn "Assembly gate: ${id} FAILED (> ${params.max_contigs} contigs) — excluded downstream" }

    ch_pass = ch_contigs
        .join(ASSEMBLY_QC.out.status)
        .filter { id, contigs, status -> status == 'PASS' }
        .map    { id, contigs, status -> tuple(id, contigs) }

    // Species per sample: Bracken top species, else --species, else 'unknown'
    if (ch_species == null)
        ch_species = ch_contigs.map { id, contigs -> tuple(id, params.species ?: 'unknown') }
    ch_pass_sp = ch_pass.join(ch_species)          // (id, contigs, organism)

    QC_SUMMARY(ch_sp_rows.collect().ifEmpty([no_file]), ch_qc_rows.collect())

    // 7. Per-genome analyses
    if (ch_ref) {
        FASTANI(ch_pass, ch_ref)
        SNIPPY(ch_pass, ch_ref)
        SNIPPY.out.dir.collect()
            .filter { it.size() >= 2 }
            .set { ch_snippy_dirs }
        SNIPPY_CORE(ch_snippy_dirs, ch_ref)
    }

    PROKKA(ch_pass_sp, params.prokka_proteins ? file(params.prokka_proteins) : no_file)
    ch_mqc = ch_mqc.mix(PROKKA.out.dir.map { it[1] })

    if (params.card_json) {
        RGI_LOAD(Channel.value(file(params.card_json)))
        RGI_MAIN(ch_pass, RGI_LOAD.out.db)
    }

    // Kaptive — only for species with a Kaptive database
    if (params.kaptive != 'off') {
        ch_kaptive_in = ch_pass_sp.flatMap { id, contigs, organism ->
            def dbs = params.kaptive == 'auto'
                ? kaptive_dbs_for(organism)
                : params.kaptive.toString().tokenize(',')*.trim()
            if (!dbs)
                log.info "Kaptive: ${id} skipped — '${organism}' is not A. baumannii or K. pneumoniae species complex"
            dbs.collect { db -> tuple(id, contigs, db) }
        }
        KAPTIVE(ch_kaptive_in)
        KAPTIVE.out.tsv
            .collectFile(keepHeader: true, sort: true, storeDir: "${params.outdir}/14_kaptive") { db, f -> ["kaptive_${db}.tsv", f] }
    }

    // 8. Pangenome (>= 2 genomes)
    if (!params.skip_pangenome) {
        PROKKA.out.gff.collect()
            .filter { gffs ->
                if (gffs.size() < 2) log.warn "Pangenome skipped: ${gffs.size()} genome(s) passed QC, at least 2 required"
                gffs.size() >= 2
            }
            .set { ch_gffs }
        ROARY(ch_gffs)
    }

    // 9. Reports
    QUAST(ch_pass.map { it[1] }.collect(), params.reference ? file(params.reference) : no_file)
    ch_mqc = ch_mqc.mix(QUAST.out.results)
    MULTIQC(ch_mqc.collect())
}

workflow.onComplete {
    log.info """
    ══════════════════════════════════════════════════
    Status    : ${workflow.success ? 'SUCCESS ✔' : 'FAILED ✘'}
    Duration  : ${workflow.duration}
    Results   : ${params.outdir}
    QC gates  : ${params.outdir}/05_qc_gates/
    ══════════════════════════════════════════════════
    """.stripIndent()
}
