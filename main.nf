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
include { FETCH_REFERENCE       } from './modules/fetch_reference'
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

// ── Curated reference genomes (assets/reference_genomes.tsv) ─────────────────

def load_reference_table(String path) {
    def f = file(path)
    if (!f.exists()) error "ERROR: --reference_list not found: ${path}"
    f.readLines().drop(1).findAll { it.trim() && !it.startsWith('#') }.collect { line ->
        def c = line.split('\t', -1)
        [ key: c[0], species: c[1],
          synonyms: (c[2] == '-' ? [] : c[2].tokenize(';')*.trim()),
          strain: c[3], assembly: c[4], is_default: c[5] == 'yes',
          who: c[6], note: c[7] ]
    }
}

REF_TABLE = load_reference_table(params.reference_list.toString())

// Default reference for a species name (Bracken name or synonym), or null
def reference_for_species(String organism) {
    def s = (organism ?: '').toLowerCase().trim()
    def hits = REF_TABLE.findAll { r -> r.species.toLowerCase() == s || r.synonyms*.toLowerCase().contains(s) }
    hits.find { it.is_default } ?: (hits ? hits[0] : null)
}

// none | file | key | auto
def reference_mode() {
    if (!params.reference)                          return 'none'
    if (params.reference == 'auto')                 return 'auto'
    if (REF_TABLE.any { it.key == params.reference }) return 'key'
    return 'file'
}

def print_reference_table() {
    def w = { s, n -> s.toString().padRight(n) }
    def out = new StringBuilder("\nCurated reference genomes (${params.reference_list})\n")
    out << w('key', 17) << w('species', 30) << w('strain', 26) << w('assembly', 18) << w('default', 8) << 'WHO 2024\n'
    out << ('-' * 108) << '\n'
    REF_TABLE.each { r ->
        out << w(r.key, 17) << w(r.species, 30) << w(r.strain, 26) << w(r.assembly, 18) << w(r.is_default ? 'yes' : '', 8) << r.who << '\n'
    }
    out << '\nUse: --reference <key>  |  --reference auto  |  --reference /path/genome.fasta\n'
    log.info out.toString()
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
    if (reference_mode() == 'file' && !file(params.reference).exists())
        error "ERROR: --reference '${params.reference}' is neither a file nor a key in ${params.reference_list} (see --list_references)"
    if (workflow.stubRun) return
    check_path('kraken_db',  params.kraken_db,  'dir')
    check_path('phix',       params.phix,       'file')
    check_path('card_json',  params.card_json,  'file')
    check_path('prokka_proteins', params.prokka_proteins, 'file')
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow {
    if (params.list_references) {
        print_reference_table()
        return
    }
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
    Reference       : ${ ['none': 'none (FastANI and Snippy skipped; see 05_qc_gates/reference_suggestion.tsv)',
                          'auto': 'auto (curated reference of the most frequent species)',
                          'key' : "${params.reference} (curated)",
                          'file': params.reference][reference_mode()] }
    CARD JSON       : ${params.card_json ?: 'none (RGI skipped)'}
    Pangenome       : ${params.skip_pangenome ? 'off' : 'Roary'}
    Kaptive         : ${params.kaptive == 'auto' ? 'auto (A. baumannii / K. pneumoniae SC samples only)' : params.kaptive}
    Profile         : ${workflow.profile}
    ══════════════════════════════════════════════════════════
    """.stripIndent()

    no_file    = file("${projectDir}/assets/NO_FILE")
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

    // Reference suggestion: most frequent species → curated reference
    ch_pass_sp.map { it[2] }.toList()
        .map { orgs ->
            def lines = ["species\tn_genomes\tsuggested_key\tstrain\tassembly\tnote"]
            orgs.countBy { it }.sort { -it.value }.each { sp, n ->
                def r = sp == 'unknown' ? null : reference_for_species(sp)
                def note = r ? "use --reference ${r.key} (or auto)" :
                           (sp == 'unknown' ? 'species unknown: set --species' : 'no curated reference: supply --reference genome.fasta')
                lines << [sp, n, r?.key ?: '-', r?.strain ?: '-', r?.assembly ?: '-', note].join('\t')
                if (r) log.info "Reference suggestion: ${n} × ${sp} → ${r.key} (${r.strain}, ${r.assembly})"
            }
            lines.join('\n') + '\n'
        }
        .collectFile(name: 'reference_suggestion.tsv', storeDir: "${params.outdir}/05_qc_gates")

    // Resolve the reference → value channel (label, species|null, fasta)
    def ref_mode = reference_mode()
    ch_refsel = null
    ch_refsel_q = null
    if (ref_mode == 'file') {
        ch_refsel = Channel.value(tuple('user', null, file(params.reference)))
    } else if (ref_mode in ['key', 'auto']) {
        if (ref_mode == 'key') {
            def r = REF_TABLE.find { it.key == params.reference }
            ch_refrow = Channel.of(tuple(r.key, r.species, r.strain, r.assembly))
        } else {
            // flatMap → queue channel that simply closes empty when no reference applies
            ch_refrow = ch_pass_sp.map { it[2] }.filter { it && it != 'unknown' }.toList()
                .flatMap { orgs ->
                    if (!orgs) { log.warn "Reference auto: no species label available (use --species) — FastANI/Snippy skipped"; return [] }
                    def top = orgs.countBy { it }.max { it.value }.key
                    def r   = reference_for_species(top)
                    if (!r) { log.warn "Reference auto: no curated reference for '${top}' — supply --reference genome.fasta; FastANI/Snippy skipped"; return [] }
                    if (orgs.unique(false).size() > 1)
                        log.warn "Reference auto: mixed species ${orgs.countBy { it }}; using ${r.key} for '${top}' — Snippy runs on '${top}' samples only"
                    [ tuple(r.key, r.species, r.strain, r.assembly) ]
                }
        }
        FETCH_REFERENCE(ch_refrow)
        ch_refsel_q = FETCH_REFERENCE.out.fasta          // queue: 0 or 1 item
        ch_refsel   = ch_refsel_q.first()
    }

    QC_SUMMARY(ch_sp_rows.collect().ifEmpty([no_file]), ch_qc_rows.collect())

    // 7. Per-genome analyses
    ch_ref = ch_refsel ? ch_refsel.map { it[2] } : null
    if (ch_ref) {
        FASTANI(ch_pass, ch_ref)
        // Core SNPs only within the reference species (unknown labels are kept)
        ch_snippy_in = ch_pass_sp.combine(ch_refsel)
            .filter { id, contigs, org, label, rsp, fasta ->
                def ok = !rsp || org == 'unknown' || org.equalsIgnoreCase(rsp)
                if (!ok) log.info "Snippy: ${id} ('${org}') excluded — reference is ${rsp}"
                ok
            }
            .map { id, contigs, org, label, rsp, fasta -> tuple(id, contigs) }
        SNIPPY(ch_snippy_in, ch_ref)
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
    ch_quast_ref = ch_refsel_q ? ch_refsel_q.map { it[2] }.concat(Channel.of(no_file)).first()
                 : (ch_ref ?: Channel.value(no_file))
    QUAST(ch_pass.map { it[1] }.collect(), ch_quast_ref)
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
