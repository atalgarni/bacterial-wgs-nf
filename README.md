# bacterial-wgs-nf

**From raw Illumina reads to annotated, quality-checked bacterial genomes, with a pangenome, a core-SNP alignment and AMR genes, in one Nextflow DSL2 workflow, for any bacterial species.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`bacterial-wgs-nf` processes paired-end short reads from bacterial **isolates**, one genome per sample. For every sample it:

1. cleans and trims the reads;
2. checks that the culture is pure and identifies its species;
3. assembles the genome and rejects poor assemblies;
4. annotates the genes and finds antimicrobial-resistance (AMR) genes.

It then compares the genomes that passed QC with each other, building a pangenome and a core-genome SNP alignment. A per-sample QC table records why any sample was excluded.

The pipeline is **species-agnostic**: nothing in it is tied to a particular genus. The organism is identified per sample from the reads and is used to label the annotation. The one exception is **Kaptive capsule/O-antigen typing**. Kaptive has databases for only two organisms, so it runs only on those (see [§3](#3-species-handling-and-the-kaptive-exception)).

---

## Contents

- [Quick start](#quick-start)
1. [What the pipeline does](#1-what-the-pipeline-does)
2. [Two ways to start: reads or assemblies](#2-two-ways-to-start-reads-or-assemblies)
3. [Species handling and the Kaptive exception](#3-species-handling-and-the-kaptive-exception)
4. [The steps in detail](#4-the-steps-in-detail)
5. [Quality-control gates](#5-quality-control-gates)
6. [Installation](#6-installation)
7. [Reference databases](#7-reference-databases)
8. [Usage recipes](#8-usage-recipes)
9. [All options](#9-all-options)
10. [Output](#10-output)
11. [Interpreting the results](#11-interpreting-the-results)
12. [Computational requirements](#12-computational-requirements)
13. [Software versions](#13-software-versions)
14. [Testing](#14-testing)
15. [Troubleshooting](#15-troubleshooting)
16. [Status and roadmap](#16-status-and-roadmap)
17. [Provenance and changes from the original scripts](#17-provenance-and-changes-from-the-original-scripts)
18. [Citation](#18-citation)
19. [References](#19-references)

---

## Quick start

```bash
# 1. Nextflow (Java 17+) in its own environment
mamba create -n nextflow -c bioconda -c conda-forge nextflow openjdk=17
mamba activate nextflow

# 2. Check the installation on the bundled test data (seconds, no tools needed)
git clone https://github.com/atalgarni/bacterial-wgs-nf.git && cd bacterial-wgs-nf
nextflow run main.nf -profile test -stub

# 3. Real run: only --reads and --kraken_db are required
nextflow run main.nf -profile conda \
    --reads 'reads/*_R{1,2}.fastq.gz' \
    --kraken_db /db/k2_standard_08gb \
    --card_json /db/card.json \
    --reference auto \
    --outdir results
```

Open `results/05_qc_gates/` first to see which samples passed, then `results/13_multiqc/multiqc_report.html`.

---

## 1. What the pipeline does

### 1.1 Flow

```
 READS ENTRY (--reads)                                     ASSEMBLY ENTRY (--assemblies)
 ───────────────────────                                   ─────────────────────────────
 R1/R2 FASTQ                                                FASTA per genome
   │                                                            │
   ├─ 1  PhiX removal ........ BWA-MEM + SAMtools   [--phix]    │
   ├─ 2  Trimming + QC ....... Trim Galore | fastp              │
   ├─ 3  Classification ...... Kraken2                          │
   ├─ 4  Abundance ........... Bracken                          │
   ├─ 5  ◆ SPECIES GATE ...... purity (+ optional target)       │
   │        └─ top species → per-sample organism label ──┐      │   organism = --species
   └─ 6  De novo assembly .... SKESA                     │      │   (or "unknown")
                    │                                    │      │
                    └──────────────┬─────────────────────┴──────┘
                                   │
                         7  ◆ ASSEMBLY GATE ........ contig count (+ N50, size, GC%)
                                   │ genomes that pass
      ┌───────────┬───────────┬────┴──────┬───────────┬────────────┬────────────┐
      ▼           ▼           ▼           ▼           ▼            ▼            ▼
  8 FastANI   9 Prokka    10 Snippy   11 RGI      12 QUAST     14 Kaptive    (per genome)
  [--reference] annotation [--reference] [--card_json]          A. baumannii /
                  │           │                                 K. pneumoniae SC only
                  ▼           ▼
           9b Roary      10b snippy-core                                     (all genomes)
           pangenome     core-SNP alignment
           (≥ 2 genomes) (≥ 2 genomes)
                                   │
                         13 MultiQC report + QC tables
```

### 1.2 Steps at a glance

| # | Step | Tool | What it gives you | Runs when | Output folder |
|---|---|---|---|---|---|
| 1 | PhiX removal | BWA-MEM, SAMtools | Reads without the Illumina PhiX spike-in | `--phix` given | `01_phix_removed` |
| 2 | Trimming + read QC | Trim Galore (default) or fastp | Adapter- and quality-trimmed reads, QC reports | always (reads entry) | `02_trim_galore` / `02_fastp` |
| 3 | Classification | Kraken2 | Read-level taxonomy report | unless `--skip_species_check` | `03_kraken2` |
| 4 | Abundance | Bracken | Species abundance per sample | as step 3 | `04_bracken` |
| 5 | **Species gate** | built-in | PASS/FAIL on purity; the top species is used as the organism label | as step 3 | `05_qc_gates` |
| 6 | Assembly | SKESA | Draft genome (contigs) | reads entry | `06_assembly` |
| 7 | **Assembly gate** | built-in | PASS/FAIL on fragmentation; N50, size, GC% | always | `05_qc_gates` |
| 8 | ANI | FastANI | Identity to the reference genome | `--reference` set (file, key or `auto`) | `07_fastani` |
| 9 | Annotation | Prokka | Genes, proteins, GFF/GenBank | always | `08_prokka` |
| 9b | Pangenome | Roary | Gene presence/absence, core-gene alignment | ≥ 2 genomes pass | `09_pangenome` |
| 10 | Variant calling | Snippy | SNPs of each genome against the reference (same-species genomes only with a curated reference) | `--reference` set | `10_snippy` |
| 10b | Core SNPs | snippy-core | Core-genome SNP alignment for phylogeny | ≥ 2 genomes + reference | `10_snippy/core` |
| 11 | AMR genes | RGI + CARD | Resistance genes and mutations | `--card_json` given | `11_rgi` |
| 12 | Assembly QC | QUAST | Assembly statistics report | always | `12_quast` |
| 13 | Report | MultiQC | One HTML page with all QC | always | `13_multiqc` |
| 14 | **K / O(C) locus typing** | Kaptive 3 | Capsule and O-antigen/outer-core locus | sample is *A. baumannii* or *K. pneumoniae* SC (auto) | `14_kaptive` |

Optional steps are switched on simply by supplying their input. Nothing breaks when an input is missing: the step is skipped and the run summary says so.

---

## 2. Two ways to start: reads or assemblies

| | `--reads` (default) | `--assemblies` |
|---|---|---|
| Input | Paired FASTQ: `'reads/*_R{1,2}.fastq.gz'` | FASTA: `'genomes/*.fasta'` |
| Steps 1–6 (cleaning, species, assembly) | run | skipped |
| Organism label | Bracken top species, per sample | `--species "Genus species"` (same for all), else unknown |
| Kraken2 database | required (unless `--skip_species_check`) | not needed |
| Use it for | New sequencing runs | Genomes from NCBI or earlier projects, or re-running comparative steps |

Sample names come from the file names. For reads, the name is the part before `_R1`/`_R2`. For assemblies, it is the file name without `.fasta`/`.fa`/`.fna` and without a trailing `_contigs`.

---

## 3. Species handling and the Kaptive exception

### 3.1 How the pipeline knows the organism

No species is hard-coded. For each sample:

1. **Reads entry:** Kraken2 + Bracken identify the most abundant species. That name (e.g. `Klebsiella pneumoniae`) becomes the sample's **organism label**.
2. **Assembly entry**, or `--skip_species_check`: the label is `--species` if you give it, otherwise `unknown`.

The label is used in two places:

| Where | How |
|---|---|
| Prokka | `--genus` / `--species` in the GenBank/GFF headers. You can override it with `--prokka_genus` / `--prokka_species`; `unknown` uses Prokka's generic defaults. Annotation itself (Prodigal + generic Prokka databases) works the same for every bacterium. |
| Kaptive (auto mode) | Decides whether Kaptive runs and with which database ([§3.2](#32-kaptive-works-for-two-organisms-only)). |

The species gate checks **purity** by default. The top species must make up ≥ 93 % of reads, whatever it is. To also require a particular species, add `--target_species "Genus species"`.

### 3.2 Kaptive works for two organisms only

> **Important:** Kaptive types the surface-polysaccharide biosynthesis loci by comparing the assembly with curated reference loci. Curated databases exist **only** for the organisms below, so for any other species Kaptive is **not run**, and its absence from the output is expected.

| Organism | Kaptive database (keyword) | Locus typed | Example result |
|---|---|---|---|
| *Acinetobacter baumannii* | `ab_k` | Capsular polysaccharide K locus (**KL**) | KL2, KL22 |
| | `ab_o` | Lipooligosaccharide outer-core locus (**OCL**) | OCL1, OCL3 |
| *Klebsiella pneumoniae* species complex: *K. pneumoniae*, *K. quasipneumoniae*, *K. variicola*, *K. quasivariicola*, *K. africana* | `kpsc_k` | Capsule K locus (**KL**) | KL64, KL47 |
| | `kpsc_o` | O-antigen (LPS) locus (**OL**) | O1/O2v1, O3b |

`--kaptive` controls the step:

| Value | Behaviour |
|---|---|
| `auto` (default) | Per sample: *A. baumannii* gets `ab_k` + `ab_o`; KpSC species get `kpsc_k` + `kpsc_o`; any other species is skipped, with a log line naming the sample and its species. |
| `off` | Never run Kaptive. |
| `ab_k,ab_o` (any comma list of the four keywords) | Run these databases on **every** sample that passes the assembly gate. Use this with `--assemblies` when no Bracken label exists and you know all genomes are *A. baumannii* (or KpSC). Do not use it for mixed-species sets: other species will come out as *Untypeable*. |

Do not over-interpret a Kaptive call from the wrong organism. The loci are species-specific, and a "match" in another genus is meaningless. This is why `auto` mode depends on the species label.

### 3.3 Choosing a reference genome: curated library and `--reference auto`

FastANI, Snippy and QUAST need a reference genome, and core-SNP analysis is only meaningful against a closed genome of the **same species**. To spare users from searching NCBI, the pipeline ships a curated table, `assets/reference_genomes.tsv`. It lists widely used, complete reference strains of the clinically most important, virulent and antibiotic-resistant bacteria: the ESKAPE organisms, the WHO Bacterial Priority Pathogens List 2024 [34], and other major pathogens. Each entry has a short **key** (e.g. `KPNE_HS11286`) that you can pass to `--reference`.

`--reference` accepts three kinds of value:

| Value | What happens |
|---|---|
| `/path/genome.fasta` | Your own reference is used for all samples (previous behaviour). |
| a key, e.g. `ABAU_ACICU` | That curated genome is downloaded from NCBI and used. |
| `auto` | The pipeline takes the **most frequent species** among the genomes that passed QC (from the species step), looks up its default key and downloads that genome. |

How a curated reference is fetched (`FETCH_REFERENCE`):

1. `datasets download genome accession <GCF_…>` (NCBI Datasets CLI [35]) downloads the assembly. This needs internet access on the machine that runs the pipeline.
2. **Safety check:** the organism name recorded by NCBI must match the expected species. If it does not, the run stops with an error, so a wrong accession can never silently become a reference.
3. Only the **chromosome(s)** are kept and plasmids are dropped, so that SNP coordinates are chromosomal. *Vibrio*, *Burkholderia* and *Brucella* keep both chromosomes.
4. The FASTA and an `.info.tsv` (accession, NCBI organism name, length) are stored in `--reference_cache` (default `./reference_cache/<key>/`) and reused by later runs without downloading again.

With a curated reference (key or `auto`), **Snippy runs only on genomes of the reference's species**. A *K. pneumoniae* genome is never aligned to an *A. baumannii* reference; the log lists the excluded samples. Genomes labelled `unknown` are kept. FastANI still runs on every genome, so off-species genomes stand out by their low ANI.

**Suggestions are always written**, even without `--reference`: `05_qc_gates/reference_suggestion.tsv` lists every species found, the number of genomes and the suggested key. You can therefore run once, look at the suggestion, and re-run with `-resume --reference <key>`.

```bash
nextflow run main.nf --list_references        # print the curated table and exit
```

Curated references (v1.0.0):

| Key | Species | Strain | Assembly | Default | WHO BPPL 2024 | Why this strain |
|---|---|---|---|---|---|---|
| `ABAU_ACICU` | *Acinetobacter baumannii* | ACICU | GCF_000018445.1 | ✔ | critical | Global clone 2 (ST2-Pasteur), carbapenem-resistant; the dominant CRAB lineage worldwide |
| `ABAU_AB5075` | *Acinetobacter baumannii* | AB5075-UW | GCF_000963815.1 |  | critical | Global clone 1, multidrug-resistant and virulent model strain |
| `KPNE_HS11286` | *Klebsiella pneumoniae* | HS11286 | GCF_000240185.1 | ✔ | critical | ST11 KPC-2 carbapenem-resistant clone; NCBI RefSeq reference genome |
| `KPNE_NTUHK2044` | *Klebsiella pneumoniae* | NTUH-K2044 | GCF_000009885.1 |  | critical | Hypervirulent K1 ST23 liver-abscess strain |
| `ECOL_MG1655` | *Escherichia coli* | K-12 substr. MG1655 | GCF_000005845.2 | ✔ | critical | Canonical E. coli reference (non-pathogenic; use for core-SNP coordinates) |
| `ECLO_ATCC13047` | *Enterobacter cloacae* | ATCC 13047 | GCF_000025565.1 | ✔ | critical | Type strain of E. cloacae; complex members are frequent carbapenemase carriers |
| `KAER_KCTC2190` | *Klebsiella aerogenes* | KCTC 2190 | GCF_000215745.1 | ✔ | critical | Complete reference of K. aerogenes (formerly Enterobacter aerogenes) |
| `SMAR_DB11` | *Serratia marcescens* | Db11 | GCF_000513215.1 | ✔ | critical | Closed reference genome of S. marcescens |
| `PMIR_HI4320` | *Proteus mirabilis* | HI4320 | GCF_000069965.1 | ✔ | critical | Uropathogenic reference strain |
| `MTUB_H37RV` | *Mycobacterium tuberculosis* | H37Rv | GCF_000195955.2 | ✔ | critical | Universal M. tuberculosis reference for variant calling and resistance-mutation coordinates |
| `PAER_PAO1` | *Pseudomonas aeruginosa* | PAO1 | GCF_000006765.1 | ✔ | high | Canonical P. aeruginosa reference |
| `PAER_PA14` | *Pseudomonas aeruginosa* | UCBPP-PA14 | GCF_000014625.1 |  | high | Highly virulent clinical strain |
| `SAUR_NCTC8325` | *Staphylococcus aureus* | NCTC 8325 | GCF_000013425.1 | ✔ | high | Canonical S. aureus reference (CC8) |
| `SAUR_USA300` | *Staphylococcus aureus* | USA300_FPR3757 | GCF_000013465.1 |  | high | Community-associated MRSA, ST8 USA300 |
| `EFAM_AUS0004` | *Enterococcus faecium* | Aus0004 | GCF_000250945.1 | ✔ | high | Vancomycin-resistant hospital-adapted clade A1 (CC17) strain |
| `EFAE_V583` | *Enterococcus faecalis* | V583 | GCF_000007785.1 | ✔ |  | First sequenced vancomycin-resistant E. faecalis |
| `SENT_LT2` | *Salmonella enterica* | serovar Typhimurium LT2 | GCF_000006945.2 | ✔ | high | Canonical non-typhoidal Salmonella reference |
| `SENT_CT18` | *Salmonella enterica* | serovar Typhi CT18 | GCF_000195995.1 |  | high | Multidrug-resistant S. Typhi reference |
| `SFLE_2A301` | *Shigella flexneri* | 2a str. 301 | GCF_000006925.2 | ✔ | high | S. flexneri serotype 2a reference |
| `SSON_SS046` | *Shigella sonnei* | Ss046 | GCF_000092525.1 | ✔ | high | S. sonnei reference |
| `NGON_FA1090` | *Neisseria gonorrhoeae* | FA 1090 | GCF_000006845.1 | ✔ | high | Canonical N. gonorrhoeae reference |
| `NMEN_MC58` | *Neisseria meningitidis* | MC58 | GCF_000008805.1 | ✔ |  | Serogroup B reference |
| `SPNE_TIGR4` | *Streptococcus pneumoniae* | TIGR4 | GCF_000006885.1 | ✔ | medium | Virulent serotype 4 reference |
| `SPYO_SF370` | *Streptococcus pyogenes* | M1 GAS SF370 | GCF_000006785.2 | ✔ | medium | emm1 group A Streptococcus reference |
| `SAGA_2603VR` | *Streptococcus agalactiae* | 2603V/R | GCF_000007265.1 | ✔ | medium | Serotype V group B Streptococcus reference |
| `HINF_RDKW20` | *Haemophilus influenzae* | Rd KW20 | GCF_000027305.1 | ✔ | medium | First sequenced bacterial genome; H. influenzae reference |
| `CJEJ_NCTC11168` | *Campylobacter jejuni* | NCTC 11168 | GCF_000009085.1 | ✔ |  | Canonical C. jejuni reference (WHO 2017 high priority) |
| `HPYL_26695` | *Helicobacter pylori* | 26695 | GCF_000008525.1 | ✔ |  | Canonical H. pylori reference (WHO 2017 high priority) |
| `CDIF_630` | *Clostridioides difficile* | 630 | GCF_000009205.2 | ✔ |  | Multidrug-resistant epidemic strain; C. difficile reference |
| `LMON_EGDE` | *Listeria monocytogenes* | EGD-e | GCF_000196035.1 | ✔ |  | Canonical L. monocytogenes reference |
| `VCHO_N16961` | *Vibrio cholerae* | O1 El Tor N16961 | GCF_000006745.1 | ✔ |  | Seventh-pandemic reference (two chromosomes) |
| `SMAL_K279A` | *Stenotrophomonas maltophilia* | K279a | GCF_000072485.1 | ✔ |  | Clinical multidrug-resistant reference |
| `BPSE_K96243` | *Burkholderia pseudomallei* | K96243 | GCF_000011545.1 | ✔ |  | Melioidosis reference (two chromosomes) |
| `BMEL_16M` | *Brucella melitensis* | 16M | GCF_000007125.1 | ✔ |  | Brucellosis reference (two chromosomes); endemic in the Middle East |
| `LPNE_PHIL1` | *Legionella pneumophila* | Philadelphia 1 | GCF_000008485.1 | ✔ |  | Original outbreak strain; L. pneumophila reference |
| `MABS_ATCC19977` | *Mycobacteroides abscessus* | ATCC 19977 | GCF_000069185.1 | ✔ |  | Type strain of the intrinsically multidrug-resistant M. abscessus |

*How the defaults were chosen:* for each species, the default is the strain most commonly used as the reference in the literature and in NCBI RefSeq, and it is a closed genome. Where a clinically dominant resistant clone differs from the classic laboratory strain, the clone was preferred. For example, *A. baumannii* ACICU (global clone 2, the dominant carbapenem-resistant lineage) and *K. pneumoniae* HS11286 (ST11, KPC-2) are the defaults. Non-default rows (e.g. `ABAU_AB5075`, `KPNE_NTUHK2044`, `SAUR_USA300`) are alternatives for other lineages; select them by key.

*Adding or changing references:* edit the TSV, or point `--reference_list` to your own table with the same columns. The best reference for an outbreak is a closed genome from the **same lineage (ST/clonal complex)** as your isolates. A distant reference within the species lowers the aligned core and can hide recombination. If your isolates belong to a lineage not represented here, choose a closed RefSeq genome of that ST and pass it as a file.

---

## 4. The steps in detail

### 4.1 PhiX removal (optional): BWA-MEM + SAMtools

Illumina runs spike in PhiX174 as a control, and leftover PhiX reads become small contaminating contigs. Reads are mapped to PhiX with BWA-MEM [3,4], and only pairs in which **both** mates are unmapped are kept [5]:

```bash
bwa mem -t $CPUS phix.fa R1.fq.gz R2.fq.gz \
  | samtools view -b -f 12 -F 256 | samtools sort -n \
  | samtools fastq -1 S_clean_R1.fastq.gz -2 S_clean_R2.fastq.gz -0 /dev/null -s /dev/null -n
```

* `-f 12` keeps a read only when it and its mate are both unmapped, so pairs stay intact.
* `-F 256` drops secondary alignments.
* `sort -n` keeps R1 and R2 in sync.

`S_phix.tsv` reports raw pairs, kept pairs and PhiX pairs removed.

### 4.2 Read trimming: Trim Galore (default) or fastp

**Trim Galore** [6] wraps Cutadapt [7]. It:

* trims low-quality 3′ ends (Phred < `--trim_quality`, default 20);
* auto-detects and removes the adapter (Illumina, Nextera or small-RNA);
* discards a pair if either read is shorter than `--trim_min_length` (default 36 bp);
* runs FastQC [8] on the trimmed reads.

```bash
trim_galore --paired --cores $N --quality 20 --length 36 --basename S --gzip --fastqc R1.fq.gz R2.fq.gz
```

`--cores` is set to the task's CPUs minus 4, limited to 1–8, because each Trim Galore core starts about three threads.

**fastp** [9] (`--trimmer fastp`) does adapter detection, quality trimming and filtering in one fast pass. The CRAB study [33] used fastp, so choose it to reproduce that analysis exactly. Both trimmers give comparable input for SKESA. Extra options go through `--trim_galore_args` or `--fastp_args`.

### 4.3 Species identification: Kraken2 + Bracken

Kraken2 [10] classifies every read against `--kraken_db`. Bracken [11] then redistributes reads that Kraken2 could only assign to genus or higher, giving species-level abundance:

```bash
kraken2 --db $DB --paired --report S.kraken2.report R1 R2 > /dev/null
bracken -d $DB -i S.kraken2.report -o S.bracken -w S_bracken.report -r $READ_LEN -l S
```

`--read_len` must match a `database{N}mers.kmer_distrib` file in the database (e.g. 100, 150).

### 4.4 De novo assembly: SKESA

SKESA [12] is a conservative de Bruijn assembler built for bacterial isolates; NCBI Pathogen Detection uses it. It does not extend through ambiguous repeats, so it produces few misassemblies, which suits SNP-based outbreak work. Contigs shorter than `--min_contig` (500 bp) are dropped.

```bash
skesa --reads R1,R2 --cores $CPUS --memory $GB --min_contig 500 > S_contigs.fasta
```

### 4.5 ANI to a reference: FastANI

With a reference (file, key or `auto`), FastANI [13] computes the whole-genome average nucleotide identity of each assembly to the reference. ≥ 95–96 % ANI means the same species. This is an assembly-based check that is independent of Kraken2.

### 4.6 Annotation: Prokka

Prokka [14] predicts:

* CDS with Prodigal [15];
* rRNA with Barrnap;
* tRNA/tmRNA with ARAGORN [16];
* ncRNA with Infernal/Rfam [17,18].

It then assigns function by a hierarchical BLAST+/HMMER search [19]:

```bash
prokka --outdir S --prefix S --strain S \
       --genus <label> --species <label> [--usegenus] [--proteins trusted.faa] \
       --addgenes --centre CDC --compliant --evalue 1e-6 --rfam --cpus $CPUS S_contigs.fasta
```

* `--proteins` (via `--prokka_proteins`) gives priority to a curated proteome of your organism. This improves the consistency of gene names, which matters for the pangenome.
* `--usegenus` (via `--prokka_usegenus`) only helps if a Prokka genus database exists for that genus, so it is off by default.

### 4.7 Pangenome: Roary

The Prokka GFF3 files of all passing genomes are clustered with Roary [20] at `--roary_identity` (95 %) BLASTP identity. A gene is *core* if it is present in ≥ `--roary_core` (99 %) of isolates.

The core genes are aligned and concatenated (`-e`). By default this uses PRANK [21], which is codon-aware and slow. `--roary_mafft` uses MAFFT [22] instead, which is much faster.

```bash
roary -f roary -e [--mafft] -i 95 -cd 99 -p $CPUS *.gff
```

Roary runs only when at least 2 genomes pass. Use genomes of **one species**: a pangenome across species is not meaningful.

### 4.8 Core-genome SNPs: Snippy + snippy-core

With a reference, each assembly is aligned to the reference with Snippy [23] in contig mode (`--ctgs`). Snippy uses BWA-MEM for mapping and FreeBayes [24] for variant calling.

`snippy-core` combines the samples into:

* `core.aln`: variable core sites, the input for phylogeny;
* `core.full.aln`: the whole-genome alignment;
* `core.txt`: per-sample alignment statistics.

Use a closed reference of the **same species** as your isolates, ideally the same lineage ([§3.3](#33-choosing-a-reference-genome-curated-library-and---reference-auto)).

### 4.9 AMR: RGI + CARD

With `--card_json`, CARD [25] is loaded locally once. Each assembly is then screened with RGI `main`, which uses DIAMOND [26] for homology and reports *Perfect* and *Strict* hits. RGI covers acquired genes and known resistance mutations for any bacterium.

### 4.10 Capsule / O-antigen typing: Kaptive 3 (species-restricted)

Kaptive 3 [27] aligns the assembly to curated reference loci [28–30] and reports:

* the best-matching locus;
* the predicted phenotype type;
* a confidence call (*Typeable* / *Untypeable*);
* problems found: missing, truncated or extra genes, or a locus split across contigs.

```bash
kaptive assembly ab_k S_contigs.fasta -o S_ab_k.tsv -t $CPUS
```

When no locus gene aligns, the sample is still listed, as `none / Untypeable`, so it is not silently missing. See [§3.2](#32-kaptive-works-for-two-organisms-only) for which samples are typed.

### 4.11 Reporting: QUAST + MultiQC

QUAST [31] summarises all passing assemblies, against the reference when one is given. MultiQC [32] combines the Trim Galore/FastQC (or fastp), Kraken2, Bracken, QUAST and Prokka reports into one HTML page.

---

## 5. Quality-control gates

### 5.1 Species gate (after Bracken)

A sample **passes** when:

1. its most abundant species has ≥ `--min_species_fraction` (default 0.93) of reads, **and**
2. if `--target_species` is set, that top species is the target (case-insensitive).

Samples that fail are **not assembled**. `05_qc_gates/species_check.tsv` explains every decision:

| sample | top_species | top_fraction | target_species | target_fraction | status | reason |
|---|---|---|---|---|---|---|
| S1 | Klebsiella pneumoniae | 0.9712 | - | - | PASS | - |
| S2 | Escherichia coli | 0.6120 | - | - | FAIL | top species below threshold |
| S3 | Klebsiella pneumoniae | 0.9500 | Acinetobacter baumannii | 0.0010 | FAIL | top species is not the target |

**Tuning for your organism.** Kraken2 cannot always separate closely related species from short reads. Examples are *E. coli* vs *Shigella*, the *Enterobacter cloacae* complex, the *Streptococcus mitis* group and the *Bacillus cereus* group. In these groups even a pure culture can score below 0.93 at species level. If pure isolates fail, you have three options:

* lower `--min_species_fraction` (e.g. 0.80);
* gate at genus level with `--bracken_level G`. The label then contains only the genus, so Kaptive `auto` will not run;
* bypass the gate with `--skip_species_check --species "Genus species"`.

### 5.2 Assembly gate (after assembly)

An assembly **passes** when it has ≤ `--max_contigs` (default 200) contigs of at least `--min_contig` bp. Heavy fragmentation usually means low coverage, a mixed culture or a failed library. Such genomes inflate the accessory genome and add false SNPs.

`05_qc_gates/assembly_qc.tsv` lists contigs, total length, largest contig, N50 and GC% for every sample. Compare total length and GC% with the expected values for your species. For genomes rich in IS elements or plasmids, a higher `--max_contigs` may be appropriate.

---

## 6. Installation

1. **Nextflow ≥ 24.04 with Java 17+:**
   ```bash
   mamba create -n nextflow -c bioconda -c conda-forge nextflow openjdk=17
   ```
2. **The pipeline:** `git clone https://github.com/atalgarni/bacterial-wgs-nf.git`, or run it directly with `nextflow run atalgarni/bacterial-wgs-nf -r v1.0.0`.
3. **The tools:** nothing to install by hand.
   * With `-profile conda`, Nextflow builds one Bioconda [2] environment per tool on first use, using mamba, and caches it under `work/conda/`. The first run therefore takes longer.
   * `-profile docker` uses Seqera Wave containers built from the same specs. *This profile is untested.*

---

## 7. Reference databases

No databases are bundled. Download them once and reuse them.

| Database | Needed for | How to get it |
|---|---|---|
| **Kraken2/Bracken** (required for the reads entry) | Species gate, organism label | A prebuilt index from <https://benlangmead.github.io/aws-indexes/k2>, e.g. *Standard-8* (8 GB RAM) or *Standard-16*. These include the Bracken `kmer_distrib` files. |
| **PhiX174** (optional) | Step 1 | NCBI RefSeq `NC_001422.1` (FASTA) |
| **Reference genome** (optional) | FastANI, Snippy, QUAST | **Nothing to download by hand** for the 31 curated species (36 genomes): use `--reference auto` or a key ([§3.3](#33-choosing-a-reference-genome-curated-library-and---reference-auto)), and it is fetched from NCBI once and cached. For other species, pass a closed genome FASTA. |
| **CARD** (optional) | RGI | `wget https://card.mcmaster.ca/latest/data && tar -xjf data ./card.json` |
| **Trusted proteins** (optional) | Prokka | A protein FASTA of a well-annotated strain of your species |
| Kaptive | Step 14 | **Nothing to download**: the four databases ship with Kaptive 3 |

All supplied paths are checked at launch, before any job starts.

---

## 8. Usage recipes

**Any species, minimum inputs.** QC, assembly, annotation and pangenome:

```bash
nextflow run main.nf -profile conda \
  --reads 'reads/*_R{1,2}.fastq.gz' --kraken_db db/k2_standard_08gb
```

**Outbreak investigation of one species.** Requires the target species and adds SNPs, AMR and Kaptive (here *K. pneumoniae*):

```bash
nextflow run main.nf -profile conda \
  --reads 'reads/*_R{1,2}.fastq.gz' --kraken_db db/k2_standard_08gb \
  --target_species "Klebsiella pneumoniae" \
  --reference KPNE_HS11286 --card_json db/card.json \
  --phix db/phix.fa --outdir kp_outbreak
```

**Carbapenem-resistant *A. baumannii*.** These are the settings of the original study [33]:

```bash
nextflow run main.nf -profile conda \
  --reads 'reads/*_R{1,2}.fastq.gz' --kraken_db db/k2_standard_08gb \
  --target_species "Acinetobacter baumannii" --trimmer fastp \
  --reference ABAU_ACICU --card_json db/card.json --phix db/phix.fa
```

**Let the pipeline pick the reference** (most frequent species → curated genome):

```bash
nextflow run main.nf -profile conda \
  --reads 'reads/*_R{1,2}.fastq.gz' --kraken_db db/k2_standard_08gb \
  --reference auto --card_json db/card.json
```

**Existing assemblies.** Typed with Kaptive explicitly, because there is no Bracken label:

```bash
nextflow run main.nf -profile conda \
  --assemblies 'genomes/*.fasta' --species "Acinetobacter baumannii" \
  --card_json db/card.json
# --species makes Kaptive auto-select ab_k + ab_o and labels Prokka;
# alternatively: --kaptive ab_k,ab_o
```

**A species without a good Kraken2 resolution (e.g. *E. coli/Shigella*):**

```bash
--min_species_fraction 0.80          # or: --bracken_level G
```

**A small machine:** `--max_cpus 4 --max_memory 12.GB`

**Faster pangenome for many genomes:** `--roary_mafft`

**Resume after a stop or a parameter change** (finished tasks are reused): add `-resume`.

---

## 9. All options

Options are given as `--name value` on the command line.

### Input / output
| Option | Default | Description |
|---|---|---|
| `--reads` | – | Glob of paired FASTQ files (quote it). Required unless `--assemblies` is used. |
| `--assemblies` | – | Glob of FASTA assemblies. Skips steps 1–6. |
| `--outdir` | `results` | Output folder. |

### Read cleaning
| Option | Default | Description |
|---|---|---|
| `--phix` | – | PhiX FASTA. Leave it out to skip PhiX removal. |
| `--trimmer` | `trim_galore` | `trim_galore` or `fastp`. |
| `--trim_quality` | `20` | Trim Galore: Phred cut-off for 3′ trimming. |
| `--trim_min_length` | `36` | Trim Galore: minimum read length after trimming. |
| `--trim_galore_fastqc` | `true` | Trim Galore: run FastQC on trimmed reads. |
| `--trim_galore_args` | `""` | Extra Trim Galore options, e.g. `"--clip_R1 5 --clip_R2 5"`. |
| `--fastp_args` | `""` | Extra fastp options (with `--trimmer fastp`). |

### Species identification and gate
| Option | Default | Description |
|---|---|---|
| `--kraken_db` | – | Kraken2/Bracken database folder. |
| `--read_len` | `100` | Bracken read length. Must match a `kmer_distrib` file in the database. |
| `--bracken_level` | `S` | Taxonomic level for Bracken and the gate (`S` species, `G` genus). |
| `--min_species_fraction` | `0.93` | Minimum fraction of reads for the top taxon. |
| `--target_species` | – | Also require that the top taxon is this one, e.g. `"Klebsiella pneumoniae"`. |
| `--skip_species_check` | `false` | Skip Kraken2, Bracken and the gate. |
| `--species` | – | Organism label when there is no Bracken result (assemblies entry or skipped gate). |

### Assembly
| Option | Default | Description |
|---|---|---|
| `--min_contig` | `500` | Minimum contig length (bp). |
| `--max_contigs` | `200` | Assembly gate: maximum number of contigs. |

### Annotation
| Option | Default | Description |
|---|---|---|
| `--prokka_genus` / `--prokka_species` | – | Override the per-sample organism label. |
| `--prokka_usegenus` | `false` | Use a Prokka genus-specific database (only if one is installed for that genus). |
| `--prokka_proteins` | – | Trusted protein FASTA used first for annotation. |
| `--prokka_centre` | `CDC` | Centre tag for GenBank-compliant contig IDs. |
| `--prokka_evalue` | `1e-6` | Similarity e-value cut-off. |

### Comparative genomics and typing
| Option | Default | Description |
|---|---|---|
| `--reference` | – | Reference for FastANI, Snippy and QUAST: a FASTA path, a curated key (e.g. `KPNE_HS11286`), or `auto`. |
| `--reference_list` | `assets/reference_genomes.tsv` | Curated reference table (replace with your own, same columns). |
| `--reference_cache` | `./reference_cache` | Where downloaded references are stored and reused. |
| `--list_references` | `false` | Print the curated table and exit. |
| `--card_json` | – | CARD `card.json` (RGI). |
| `--kaptive` | `auto` | `auto`, `off`, or a comma list of `ab_k,ab_o,kpsc_k,kpsc_o` (see [§3.2](#32-kaptive-works-for-two-organisms-only)). |
| `--skip_pangenome` | `false` | Do not run Roary. |
| `--roary_identity` | `95` | Minimum BLASTP % identity. |
| `--roary_core` | `99` | Percentage of isolates a gene must be in to count as core. |
| `--roary_core_alignment` | `true` | Build the core-gene alignment. |
| `--roary_mafft` | `false` | Use MAFFT instead of PRANK (much faster). |

### Resources
| Option | Default | Description |
|---|---|---|
| `--max_cpus` | the smaller of 16 and the available cores | Upper limit per task. |
| `--max_memory` | `16.GB` | Upper limit per task. |
| `--max_time` | `48.h` | Upper limit per task. |

Nextflow's own options take a single dash: `-profile`, `-resume`, `-stub`, `-work-dir`.

---

## 10. Output

```
results/
├── 01_phix_removed/   S_phix.tsv                          PhiX read counts
├── 02_trim_galore/    *_trimming_report.txt, S_val_{1,2}_fastqc.{html,zip}
│   (02_fastp/         S_fastp.{html,json}                 with --trimmer fastp)
├── 03_kraken2/        S.kraken2.report
├── 04_bracken/        S.bracken, S_bracken.report
├── 05_qc_gates/       species_check.tsv, assembly_qc.tsv,  ← start here
│                      reference_suggestion.tsv
├── 06_assembly/       S_contigs.fasta
├── 07_fastani/        S.ani.txt
├── 08_prokka/S/       S.gff  S.gbk  S.faa  S.ffn  S.fna  S.tsv  S.txt …
├── 09_pangenome/roary/
│                      gene_presence_absence.{csv,Rtab}, summary_statistics.txt,
│                      core_gene_alignment.aln, accessory_binary_genes.fa.newick …
├── 10_snippy/S/       snps.vcf, snps.tab …
│   └── core/          core.aln, core.full.aln, core.tab, core.txt, core.vcf
├── 11_rgi/            S.txt, S.json
├── 12_quast/quast/    report.html, report.tsv …
├── 13_multiqc/        multiqc_report.html
├── 14_kaptive/        kaptive_ab_k.tsv, kaptive_ab_o.tsv, kaptive_kpsc_k.tsv, kaptive_kpsc_o.tsv
│   └── per_sample/    S_<db>.tsv
└── pipeline_info/     report.html, timeline.html, trace.txt, dag.html

reference_cache/<KEY>/ <KEY>.fasta, <KEY>.info.tsv      curated references (outside results/, reused)
```

The `14_kaptive/` summary tables exist only for databases that were actually used.

---

## 11. Interpreting the results

| File | What to look at |
|---|---|
| `05_qc_gates/*.tsv` | Every excluded sample, with the reason. Check `top_species` for unexpected organisms, and total length and GC% against your species. |
| `07_fastani/*.ani.txt` | ≥ 95 % means the same species as the reference. 80–95 % suggests a different species of the same genus, which can mean a misidentification in culture. |
| `09_pangenome/.../summary_statistics.txt` | Core / soft-core / shell / cloud gene counts. A cloud that keeps growing as genomes are added means an open pangenome. `gene_presence_absence.csv` can be loaded into Phandango or Scoary. |
| `10_snippy/core/core.txt` | Aligned fraction per sample. A low value means the reference is too distant. Use `snp-dists core.aln` for pairwise SNP distances and IQ-TREE for a tree. |
| `11_rgi/*.txt` | Prefer *Strict*/*Perfect* hits. Remember that intrinsic, chromosomal genes appear in every isolate of some species (e.g. *bla*<sub>OXA-51-like</sub> and *ampC*/ADC in *A. baumannii*, *bla*<sub>SHV</sub> in *K. pneumoniae*). Whether these are expressed (e.g. through upstream IS*Aba1*) is not reported by RGI. |
| `14_kaptive/*.tsv` | `Best match locus` (e.g. KL2, OCL1, KL64) and `Match confidence`. *Untypeable* with "problems" usually means the locus is split across contigs or interrupted by an IS element. Check `Problems`, `Coverage` and the missing/extra gene columns before reporting a novel locus. |

---

## 12. Computational requirements

These are typical per-sample figures for a genome of about 5 Mb at about 100× coverage with PE150 reads.

| Step | Resources | Time |
|---|---|---|
| Kraken2 | RAM ≈ database size (8 GB for Standard-8) | minutes |
| Trim Galore | 2–8 cores | 5–10 min (fastp: ~1–2 min) |
| SKESA | 4–8 GB | 5–15 min on 8 cores |
| Prokka | 8 cores | 5–10 min |
| RGI, Kaptive | low | 1–5 min each |
| Roary + PRANK | the bottleneck | hours for > 100 genomes; use `--roary_mafft` |

The defaults fit a 16 GB workstation. For a larger Kraken2 index (e.g. PlusPF, about 70 GB), raise `--max_memory` and the `KRAKEN2` memory in `nextflow.config`.

---

## 13. Software versions

These versions are pinned in `nextflow.config` (Bioconda).

| Tool | Version | Step |
|---|---|---|
| BWA | 0.7.18 | PhiX removal |
| SAMtools | 1.20 | PhiX removal |
| Trim Galore | 0.6.10 | trimming (default) |
| fastp | 0.23.4 | trimming (`--trimmer fastp`) |
| Kraken2 | 2.1.3 | classification |
| Bracken | 2.9 | abundance |
| SKESA | 2.5.1 | assembly |
| FastANI | 1.34 | ANI |
| Prokka | 1.14.6 | annotation |
| Roary | 3.13.0 | pangenome |
| Snippy | 4.6.0 | core SNPs |
| RGI | 6.0.3 | AMR |
| Kaptive | 3.1.0 | K / O(C) locus typing |
| QUAST | 5.2.0 | assembly QC |
| MultiQC | 1.25 | report |
| NCBI Datasets CLI | latest (conda-forge, unpinned) | reference download only |

For publications, also report the CARD version (shown by `rgi database --version --local`) and the Kraken2 database name and date.

---

## 14. Testing

```bash
nextflow run main.nf -profile test -stub
```

This uses the tiny dataset in `test_data/` and runs every branch of the workflow in a few seconds, without the real tools. The stub assigns the three test isolates to *A. baumannii*, *K. pneumoniae* and *E. coli*. You should therefore see Kaptive run with `ab_*` and `kpsc_*` databases on the first two and skip the third, with a log message. Add `--reference auto` to exercise the reference selection: the stub picks a curated reference, and Snippy runs only on the matching species.

---

## 15. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `Process requirement exceeds available memory` | Set `--max_memory` to your RAM (e.g. `12.GB`). |
| Bracken: `kmer distribution file not found` | `--read_len` has no matching `database*mers.kmer_distrib` in the database. |
| Pure isolates FAIL the species gate | The species group is hard for Kraken2 ([§5.1](#51-species-gate-after-bracken)). Lower `--min_species_fraction` or use `--bracken_level G`. |
| All samples FAIL with "top species is not the target" | The spelling in `--target_species` must match the Bracken name. Check a `.bracken` file. |
| `14_kaptive` empty or missing | No sample was *A. baumannii* or KpSC, which is expected for other species. With `--assemblies`, add `--species` or `--kaptive ab_k,ab_o`. |
| Kaptive: `none / Untypeable` | Wrong species, or the locus is fragmented in the assembly. |
| `Pangenome skipped` | Fewer than 2 genomes passed both gates. |
| No FastANI/Snippy/RGI output | `--reference` / `--card_json` was not given, or `--reference auto` found no curated reference (see the warning and `reference_suggestion.tsv`). |
| `FETCH_REFERENCE` fails to download | The machine has no internet access, or NCBI is unreachable. Download the genome elsewhere and pass it as a file. |
| `FETCH_REFERENCE`: "is '…', expected '…'" | The accession in the table does not belong to that species. Correct the row in `reference_genomes.tsv` (or report it). |
| Snippy skipped some samples | They are a different species from the curated reference (listed in the log). This is intended. |
| `--phix null` still runs PhiX removal | Leave the option out entirely instead. |
| Conda environment build hangs | Install mamba. Delete `work/conda/` and retry. |

---

## 16. Status and roadmap

**v1.0.0:**

* All branches are verified with stub runs.
* Both QC gates passed real-mode tests on synthetic data.
* The Trim Galore and Kaptive commands passed real-mode tests on the test data.
* A full end-to-end run on real isolates in this consolidated form is pending.
* The Docker profile is untested.
* The curated reference accessions and the NCBI download step have not yet been tested against live NCBI. The species check in `FETCH_REFERENCE` stops the run if an accession is wrong; please report any such row.

**Planned for v1.1:**

* MLST (PubMLST schemes, auto-detected)
* CheckM2 completeness/contamination as an additional assembly gate
* Panaroo as an alternative to Roary
* IQ-TREE phylogeny and `snp-dists` on the core alignment
* Optional Gubbins recombination masking
* Samplesheet (CSV) input

---

## 17. Provenance and changes from the original scripts

The workflow brings together steps that were first run as separate bash scripts and partial Nextflow files: PhiX removal, fastp and Kraken2 in bash; Bracken onwards in Nextflow. The carbapenem-resistant *A. baumannii* isolates from Makkah and Riyadh were analysed with those steps [33].

v1.0.0 makes them a single pipeline. It also makes the defaults species-agnostic: previously the target species and the Prokka genus were fixed to *A. baumannii*. It also fixes the following defects:

| Area | Original behaviour | v1.0.0 |
|---|---|---|
| Species gate | Read **row 2** of the unsorted Bracken table and never checked the species name, so a contaminated or wrong-species sample could pass. | Uses the most abundant species. The name is checked only when `--target_species` is set. |
| PhiX removal | `samtools view -f 4` kept orphan mates and broke R1/R2 pairing. | `-f 12 -F 256` + `sort -n`. |
| Pangenome guard | `[ $n -lt 2]` was a bash syntax error. | Genome count checked in Groovy. |
| Core SNPs | `snippy-multi` with absolute `work/` paths broke on `-resume`. | Per-sample `snippy --ctgs` + `snippy-core`. |
| Structure | Several partial `main*.nf` files with hard-coded paths. | One `main.nf`, 16 modules, validated parameters, profiles and a stub test. |

Because of these fixes, results can differ slightly from the original runs.

---

## 18. Citation

Please cite the pipeline (see `CITATION.cff`; the Zenodo DOI will be added at release) **and** the tools you used from [§19](#19-references).

> Algarni AT. *bacterial-wgs-nf: a species-agnostic Nextflow workflow for short-read bacterial genome assembly, annotation, pangenome, core-SNP and AMR analysis.* https://github.com/atalgarni/bacterial-wgs-nf

---

## 19. References

1. Di Tommaso P, Chatzou M, Floden EW, *et al.* Nextflow enables reproducible computational workflows. *Nat Biotechnol.* 2017;35:316–319. doi:10.1038/nbt.3820
2. Grüning B, Dale R, Sjödin A, *et al.* Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nat Methods.* 2018;15:475–476. doi:10.1038/s41592-018-0046-7
3. Li H, Durbin R. Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics.* 2009;25(14):1754–1760. doi:10.1093/bioinformatics/btp324
4. Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv.* 2013;1303.3997.
5. Danecek P, Bonfield JK, Liddle J, *et al.* Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. doi:10.1093/gigascience/giab008
6. Krueger F, James F, Ewels P, *et al.* FelixKrueger/TrimGalore. Zenodo. doi:10.5281/zenodo.5127898
7. Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal.* 2011;17(1):10–12. doi:10.14806/ej.17.1.200
8. Andrews S. FastQC: a quality control tool for high throughput sequence data. 2010. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
9. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics.* 2018;34(17):i884–i890. doi:10.1093/bioinformatics/bty560
10. Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biol.* 2019;20:257. doi:10.1186/s13059-019-1891-0
11. Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. *PeerJ Comput Sci.* 2017;3:e104. doi:10.7717/peerj-cs.104
12. Souvorov A, Agarwala R, Lipman DJ. SKESA: strategic k-mer extension for scrupulous assemblies. *Genome Biol.* 2018;19:153. doi:10.1186/s13059-018-1540-z
13. Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nat Commun.* 2018;9:5114. doi:10.1038/s41467-018-07641-9
14. Seemann T. Prokka: rapid prokaryotic genome annotation. *Bioinformatics.* 2014;30(14):2068–2069. doi:10.1093/bioinformatics/btu153
15. Hyatt D, Chen GL, LoCascio PF, *et al.* Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics.* 2010;11:119. doi:10.1186/1471-2105-11-119
16. Laslett D, Canback B. ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. *Nucleic Acids Res.* 2004;32(1):11–16. doi:10.1093/nar/gkh152
17. Nawrocki EP, Eddy SR. Infernal 1.1: 100-fold faster RNA homology searches. *Bioinformatics.* 2013;29(22):2933–2935. doi:10.1093/bioinformatics/btt509
18. Kalvari I, Nawrocki EP, Ontiveros-Palacios N, *et al.* Rfam 14: expanded coverage of metagenomic, viral and microRNA families. *Nucleic Acids Res.* 2021;49(D1):D192–D200. doi:10.1093/nar/gkaa1047
19. Camacho C, Coulouris G, Avagyan V, *et al.* BLAST+: architecture and applications. *BMC Bioinformatics.* 2009;10:421. doi:10.1186/1471-2105-10-421
20. Page AJ, Cummins CA, Hunt M, *et al.* Roary: rapid large-scale prokaryote pan genome analysis. *Bioinformatics.* 2015;31(22):3691–3693. doi:10.1093/bioinformatics/btv421
21. Löytynoja A. Phylogeny-aware alignment with PRANK. *Methods Mol Biol.* 2014;1079:155–170. doi:10.1007/978-1-62703-646-7_10
22. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol Biol Evol.* 2013;30(4):772–780. doi:10.1093/molbev/mst010
23. Seemann T. Snippy: fast bacterial variant calling from NGS reads. 2015. https://github.com/tseemann/snippy
24. Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. *arXiv.* 2012;1207.3907.
25. Alcock BP, Huynh W, Chalil R, *et al.* CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res.* 2023;51(D1):D690–D699. doi:10.1093/nar/gkac920
26. Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. *Nat Methods.* 2015;12:59–60. doi:10.1038/nmeth.3176
27. Stanton TD, Hetland MAK, Löhr IH, Holt KE, Wyres KL. Fast and accurate in silico antigen typing with Kaptive 3. *Microb Genom.* 2025;11(6):001428. doi:10.1099/mgen.0.001428
28. Lam MMC, Wick RR, Judd LM, Holt KE, Wyres KL. Kaptive 2.0: updated capsule and lipopolysaccharide locus typing for the *Klebsiella pneumoniae* species complex. *Microb Genom.* 2022;8(3):000800. doi:10.1099/mgen.0.000800
29. Wyres KL, Cahill SM, Holt KE, Hall RM, Kenyon JJ. Identification of *Acinetobacter baumannii* loci for capsular polysaccharide (KL) and lipooligosaccharide outer core (OCL) synthesis in genome assemblies using curated reference databases compatible with Kaptive. *Microb Genom.* 2020;6(3):e000339. doi:10.1099/mgen.0.000339
30. Cahill SM, Hall RM, Kenyon JJ. An update to the database for *Acinetobacter baumannii* capsular polysaccharide locus typing extends the extensive and diverse repertoire of genes found at and outside the K locus. *Microb Genom.* 2022;8(10):mgen000878. doi:10.1099/mgen.0.000878
31. Gurevich A, Saveliev V, Vyahhi N, Tesler G. QUAST: quality assessment tool for genome assemblies. *Bioinformatics.* 2013;29(8):1072–1075. doi:10.1093/bioinformatics/btt086
32. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics.* 2016;32(19):3047–3048. doi:10.1093/bioinformatics/btw354
33. Alquait AA, …, Algarni AT, *et al.* Phenotypic, PCR, and whole-genome sequencing characterization of carbapenem-resistant *Acinetobacter baumannii* in Makkah and Riyadh, Saudi Arabia: AST discordance and local genomic micro-clusters. *Front Microbiol.* 2026;17:1845444.
34. World Health Organization. *WHO Bacterial Priority Pathogens List, 2024: bacterial pathogens of public health importance to guide research, development and strategies to prevent and control antimicrobial resistance.* Geneva: WHO; 2024. https://www.who.int/publications/i/item/9789240093461
35. O'Leary NA, Cox E, Holmes JB, *et al.* Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. *Sci Data.* 2024;11:732. doi:10.1038/s41597-024-03571-y

---

**Author:** Abdulmalek T. Algarni ([ORCID 0000-0002-2872-3724](https://orcid.org/0000-0002-2872-3724)). Released under the MIT License.
