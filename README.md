# bacterial-wgs-nf

**Short-read bacterial whole-genome sequencing: QC, species verification, de novo assembly, annotation, pangenome, core-SNP and AMR profiling in a single Nextflow DSL2 workflow.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`bacterial-wgs-nf` takes paired-end Illumina reads (or pre-built assemblies) from bacterial isolates and produces:

* quality-controlled draft genomes
* functional annotation
* a Roary pangenome with a core-gene alignment
* a Snippy core-SNP alignment against a reference
* CARD-based antimicrobial-resistance (AMR) calls

Two explicit QC gates, one on taxonomic purity and one on assembly fragmentation, stop contaminated or poor genomes from entering the comparative analyses.

The workflow was developed for *Acinetobacter baumannii* surveillance, and its defaults reflect that. Every organism-specific setting is a parameter (target species, Prokka genus/species, reference genome), so it can be used for any bacterial species.

---

## Contents

1. [Background and provenance](#1-background-and-provenance)
2. [Workflow overview](#2-workflow-overview)
3. [Methods, step by step](#3-methods-step-by-step)
4. [Quality-control gates](#4-quality-control-gates)
5. [Installation](#5-installation)
6. [Reference databases](#6-reference-databases)
7. [Input](#7-input)
8. [Usage](#8-usage)
9. [Parameters](#9-parameters)
10. [Output](#10-output)
11. [Interpreting the results](#11-interpreting-the-results)
12. [Computational requirements](#12-computational-requirements)
13. [Software versions](#13-software-versions)
14. [Testing](#14-testing)
15. [Troubleshooting](#15-troubleshooting)
16. [Known limitations and roadmap](#16-known-limitations-and-roadmap)
17. [Citation](#17-citation)
18. [References](#18-references)

---

## 1. Background and provenance

This workflow brings together analysis steps that were first run as separate pieces:

* PhiX removal, fastp and Kraken2 as bash scripts
* Bracken onwards as partial Nextflow files

Those same steps were used for the genomic analysis of carbapenem-resistant *A. baumannii* (CRAB) isolates from Makkah and Riyadh [1].

Version 1.0.0 consolidates them into one reproducible pipeline. The order of steps is unchanged. The one tool change is read trimming: the default is now Trim Galore, and fastp, which the published analysis used, is still available with `--trimmer fastp`. Several defects in the original scripts were fixed on the way (see [§4.3](#43-fixes-relative-to-the-original-scripts)). Because of these fixes, results can differ slightly from the original bash runs. Most notably, samples that previously passed the species gate by chance can now fail it.

## 2. Workflow overview

```
                    paired-end reads (R1/R2)                pre-built assemblies
                              │                                    │ (--assemblies)
                 ┌────────────▼────────────┐                       │
  (optional) ──► │ 1  PhiX removal         │ BWA-MEM + SAMtools    │
                 └────────────┬────────────┘                       │
                 ┌────────────▼────────────┐                       │
                 │ 2  Read trimming / QC   │ Trim Galore | fastp   │
                 └────────────┬────────────┘                       │
                 ┌────────────▼────────────┐                       │
                 │ 3  Taxonomic profiling  │ Kraken2               │
                 │ 4  Abundance re-estim.  │ Bracken               │
                 │ 5  SPECIES GATE ◆       │ top species = target  │
                 └────────────┬────────────┘   and ≥ min fraction  │
                              │ PASS                               │
                 ┌────────────▼────────────┐                       │
                 │ 6  De novo assembly     │ SKESA                 │
                 └────────────┬────────────┘                       │
                              └──────────────┬─────────────────────┘
                              ┌──────────────▼──────────────┐
                              │ 7  ASSEMBLY GATE ◆          │ contigs ≤ max_contigs
                              └──────────────┬──────────────┘ (+ N50, size, GC%)
                                             │ PASS
        ┌──────────────┬──────────────┬──────┴───────┬──────────────┐
        ▼              ▼              ▼              ▼              ▼
   8 FastANI      9 Prokka       10 Snippy       11 RGI         12 QUAST
   (vs ref)       annotation     (--ctgs, ref)   CARD AMR       assembly QC
                     │              │
                     ▼              ▼
               9b Roary         10b snippy-core
               pangenome        core-SNP alignment
               (≥2 genomes)     (≥2 genomes)
                                             │
                              ┌──────────────▼──────────────┐
                              │ 13 MultiQC + QC summary     │
                              └─────────────────────────────┘
```

Some steps depend on optional inputs and are skipped cleanly when the input is missing:

| Steps | Need |
|---|---|
| 8 (FastANI) and 10 (Snippy) | `--reference` |
| 11 (RGI) | `--card_json` |

## 3. Methods, step by step

### 3.1 PhiX removal (optional): BWA-MEM + SAMtools

Illumina runs spike in PhiX174 as a sequencing control. Any PhiX reads that remain end up as small contaminating contigs.

When `--phix` is given, the reference is indexed and reads are mapped with BWA-MEM [2,3]. Only pairs in which **both** mates are unmapped are kept:

```bash
bwa index phix.fa
bwa mem -t $CPUS phix.fa R1.fq.gz R2.fq.gz \
  | samtools view -b -f 12 -F 256 \
  | samtools sort -n \
  | samtools fastq -1 clean_R1.fastq.gz -2 clean_R2.fastq.gz -0 /dev/null -s /dev/null -n
```

* `-f 12` keeps reads where the read is unmapped **and** its mate is unmapped, so only whole pairs survive. The original script used `-f 4`, which kept orphan mates and broke pairing.
* `-F 256` drops secondary alignments.
* `sort -n` puts mates next to each other, so `samtools fastq` writes synchronised R1/R2 files [4].

A per-sample `<S>_phix.tsv` records raw pairs, kept pairs and PhiX pairs removed.

### 3.2 Read trimming and QC: Trim Galore (default) or fastp

The trimmer is chosen with `--trimmer`.

**Trim Galore (default, `--trimmer trim_galore`).** Trim Galore [27] is a wrapper around Cutadapt [28]. It works in this order:

1. Trims low-quality 3′ ends (Phred < `--trim_quality`, default 20) with the BWA trimming algorithm.
2. Auto-detects the adapter (Illumina universal, Nextera or small RNA) from the first million reads and removes it. Even 1 bp of overlap at the 3′ end is removed.
3. Discards pairs in which either read is shorter than `--trim_min_length` (default 36 bp), so R1 and R2 stay synchronised.
4. Runs FastQC [29] on the trimmed reads (`--trim_galore_fastqc`, default on).

```bash
trim_galore --paired --cores $N --quality 20 --length 36 \
            --basename S --gzip --fastqc R1.fq.gz R2.fq.gz
# → S_val_1.fq.gz, S_val_2.fq.gz, *_trimming_report.txt, *_fastqc.{html,zip}
```

`--cores` is set to `task.cpus − 4`, limited to 1–8, because each Trim Galore core starts about three threads (Cutadapt plus pigz). Extra options go through `--trim_galore_args`, for example `--trim_galore_args "--clip_R1 5 --clip_R2 5"`.

**fastp (`--trimmer fastp`).** fastp [5] does adapter auto-detection for paired-end reads, sliding-window quality trimming and length filtering in a single pass. It is faster, and its HTML/JSON reports feed MultiQC. fastp's own defaults are used, and extra arguments go through `--fastp_args`. Choose fastp to reproduce the trimming exactly as in the CRAB study [1].

**Which to choose.** Both trimmers produce comparable reads for de novo assembly of bacterial isolates, and the downstream steps accept either.

* Trim Galore uses well-established Cutadapt adapter matching, trims conservatively and includes FastQC.
* fastp is several times faster and also corrects overlapping bases in paired reads.

### 3.3 Taxonomic classification: Kraken2

Trimmed reads are classified against a Kraken2 database [6] (`--kraken_db`), producing a standard report. Per-read output is discarded to save disk.

### 3.4 Abundance re-estimation: Bracken

Kraken2 assigns many reads to higher taxonomic levels when related species share k-mers. This is common in the *A. calcoaceticus–baumannii* complex. Bracken [7] uses a Bayesian model to redistribute those reads and estimate species-level abundance:

```bash
bracken -d $DB -i S.kraken2.report -o S.bracken -w S_bracken.report -r $READ_LEN -l S
```

`--read_len` must match a `databaseXmers.kmer_distrib` file inside the Kraken2 database (e.g. 100 or 150).

### 3.5 Species gate

See [§4.1](#41-species-gate).

### 3.6 De novo assembly: SKESA

SKESA [8] is a conservative de Bruijn graph assembler designed for bacterial isolates. It favours contig correctness over length and does not extend through ambiguous repeats. That makes it well suited to SNP-based outbreak work, and it is the assembler NCBI uses for Pathogen Detection.

```bash
skesa --reads R1.fq.gz,R2.fq.gz --cores $CPUS --memory $GB --min_contig 500 > S_contigs.fasta
```

### 3.7 Assembly gate

See [§4.2](#42-assembly-gate). The pipeline also computes, for each assembly:

* contig count
* total length
* largest contig
* N50
* GC%

### 3.8 Average nucleotide identity: FastANI

When a `--reference` genome is given, FastANI [9] calculates whole-genome ANI between each assembly and the reference. This confirms the Bracken call independently, from the assembly rather than the reads. The usual species boundary is an ANI of at least 95–96 %.

### 3.9 Genome annotation: Prokka

Prokka [10] predicts:

* CDS with Prodigal [11]
* rRNA with Barrnap
* tRNA/tmRNA with Aragorn [12]
* ncRNA with Infernal/Rfam [13,14], because `--rfam` is set

CDS are annotated by a hierarchical BLAST+/HMMER search [15]. The pipeline runs:

```bash
prokka --outdir $S --prefix $S --locustag $S --cpus $CPUS \
       --genus Acinetobacter --species baumannii --usegenus \
       --centre CDC --compliant --evalue 1e-6 --rfam \
       [--proteins trusted.faa] S_contigs.fasta
```

* `--usegenus` turns on the genus-specific BLAST database.
* `--proteins` (optional, via `--prokka_proteins`) gives priority to a curated protein set, such as a well-annotated reference proteome. This makes gene names more consistent across isolates.
* `--compliant --centre` produces contig IDs that are compatible with GenBank submission.

### 3.10 Pangenome: Roary

The GFF3 files from all genomes that passed the assembly gate are clustered with Roary [16]:

* Clustering uses `--roary_identity` (default 95 %) BLASTP identity.
* Core genes are those present in at least `--roary_core` (default 99 %) of isolates.

```bash
roary -f roary -e [--mafft] -i 95 -cd 99 -p $CPUS *.gff
```

With `-e`, each core gene is aligned and the alignments are concatenated into `core_gene_alignment.aln`. The aligner is one of:

* PRANK [17], the default: codon-aware but slow
* MAFFT [18], with `--roary_mafft`: much faster

Roary runs only when **2 or more** genomes are available; otherwise the step is skipped with a warning.

### 3.11 Core-SNP alignment: Snippy

When `--reference` is given, each assembly is aligned to the reference with Snippy [19] in contig mode (`--ctgs`). Snippy shreds the contigs into synthetic reads, maps them with BWA-MEM and calls variants with FreeBayes [20].

`snippy-core` then combines all samples and writes:

* `core.aln`: the core-genome SNP alignment
* `core.full.aln`: the full alignment
* `core.txt`: a per-sample SNP summary

The alignment is the input for phylogenetics (e.g. IQ-TREE or RAxML-NG) and for clustering by pairwise SNP distance (e.g. `snp-dists`).

Each sample runs as a separate task, and `snippy-core` reads directories staged into its own work folder. This replaced the original `snippy-multi` call, which wrote absolute work-directory paths into a tab file and broke whenever `-resume` changed hashes.

### 3.12 AMR profiling: RGI / CARD

When `--card_json` is supplied:

1. The CARD database [21] is loaded locally once: `rgi load --card_json card.json --local`.
2. Each assembly is screened with RGI `main` in contig mode.

RGI uses DIAMOND [22] for protein homology and reports Perfect and Strict hits. Each sample produces `<S>.txt` (tabular) and `<S>.json`.

### 3.13 Assembly QC and reporting: QUAST + MultiQC

QUAST [23] summarises every assembly that passed, against the reference when one is given. MultiQC [24] combines the Trim Galore/FastQC (or fastp), Kraken2, QUAST and Prokka reports into one HTML page. A separate `05_qc_gates/` folder holds the pipeline's own gate decisions as TSV tables.

## 4. Quality-control gates

### 4.1 Species gate

After Bracken, a sample gets **PASS** only if both conditions hold:

1. The species with the **highest** `fraction_total_reads` reaches at least `--min_species_fraction` (default **0.93**).
2. If `--target_species` is set (default *Acinetobacter baumannii*), that top species **is** the target. The match is exact but ignores case.

Samples that fail are **not assembled**. A warning is printed, and the reason is recorded in `05_qc_gates/species_check.tsv`:

| sample | top_species | top_fraction | target | target_fraction | status | reason |
|---|---|---|---|---|---|---|
| S1 | Acinetobacter baumannii | 0.9712 | Acinetobacter baumannii | 0.9712 | PASS | ok |
| S2 | Klebsiella pneumoniae | 0.9500 | Acinetobacter baumannii | 0.0210 | FAIL | top species is not the target |
| S3 | Acinetobacter baumannii | 0.8000 | Acinetobacter baumannii | 0.8000 | FAIL | top fraction below threshold |

For species other than *A. baumannii*:

* `--target_species "Genus species"`: set the target. The name must match how Bracken writes it.
* `--target_species ""`: check purity only.
* `--skip_species_check`: bypass the gate entirely. No Kraken2 database is then needed.

### 4.2 Assembly gate

An assembly gets **FAIL** when it has more than `--max_contigs` (default **200**) contigs of at least `--min_contig` bp. Failed assemblies are excluded from FastANI, Prokka, Roary, Snippy, RGI and QUAST.

Heavy fragmentation usually means low coverage, mixed cultures or failed libraries. Such genomes inflate the accessory genome and add false SNPs to the core alignment.

`05_qc_gates/assembly_qc.tsv` records the following for every sample, so borderline genomes can be reviewed by hand:

* contig count
* total length
* largest contig
* N50
* GC%

### 4.3 Fixes relative to the original scripts

| Area | Original behaviour | v1.0.0 |
|---|---|---|
| Species gate | Read **row 2** of the Bracken output, which is not sorted by abundance, and never checked the species name. A sample dominated by another organism could pass. | Takes the maximum of `fraction_total_reads` and checks the name against `--target_species`. |
| PhiX removal | `samtools view -f 4` kept singletons and broke R1/R2 synchronisation. | `-f 12 -F 256` + `sort -n`: only intact pairs are kept. |
| Pangenome guard | `[ $n -lt 2]` was a bash syntax error, so the guard never worked. | The genome count is checked in Groovy before Roary is launched. |
| Core SNPs | `snippy-multi` with absolute `work/` paths broke on `-resume` and on other machines. | Per-sample `snippy --ctgs` + `snippy-core` on staged inputs. |
| Structure | Several partial `main*.nf` files with hard-coded local paths and undefined processes; the first half lived in bash. | One `main.nf` and 14 modules, two entry points (reads / assemblies), validated parameters, conda/Docker profiles and a stub test. |

## 5. Installation

### 5.1 Nextflow

Nextflow 24.04 or later needs Java 17 or later. The simplest way to get both is a dedicated environment:

```bash
mamba create -n nextflow -c bioconda -c conda-forge nextflow openjdk=17
mamba activate nextflow
nextflow -version
```

### 5.2 Pipeline

```bash
git clone https://github.com/atalgarni/bacterial-wgs-nf.git
cd bacterial-wgs-nf
```

Or let Nextflow fetch it: `nextflow run atalgarni/bacterial-wgs-nf -r v1.0.0 ...`

### 5.3 Software

With `-profile conda`, Nextflow builds one environment per tool from Bioconda [25] on first use and caches them in `work/conda/`. The profile sets `conda.useMamba = true`, so mamba must be installed.

`-profile docker` uses Seqera Wave containers built from the same Bioconda specs. **This profile has not been tested yet.**

## 6. Reference databases

No databases are bundled. Download them once and pass their paths.

| Database | Needed for | How to obtain |
|---|---|---|
| **Kraken2 + Bracken** | species gate | A prebuilt index (e.g. *Standard-8* or *PlusPF-16*) from <https://benlangmead.github.io/aws-indexes/k2>. These already include the `database{100,150}mers.kmer_distrib` files. Extract it and pass the folder with `--kraken_db`. |
| **PhiX174** | PhiX removal | NCBI RefSeq `NC_001422.1` as FASTA, passed with `--phix`. |
| **Reference genome** | FastANI, Snippy, QUAST | A closed genome of the target species (e.g. *A. baumannii* ATCC 17978, or a locally relevant ST2 strain), passed with `--reference`. |
| **CARD** | RGI | `wget https://card.mcmaster.ca/latest/data && tar -xjf data ./card.json`, then `--card_json card.json`. |
| **Trusted proteins** (optional) | Prokka | A FASTA proteome of a well-annotated strain, passed with `--prokka_proteins`. |

All supplied paths are checked at launch. A missing path stops the run before any task starts.

## 7. Input

### 7.1 Reads (default entry)

Pass a glob that matches paired FASTQ files. The sample ID is the part before `_R1`/`_R2`:

```
reads/
├── AB001_R1.fastq.gz
├── AB001_R2.fastq.gz
├── AB002_R1.fastq.gz
└── AB002_R2.fastq.gz
```

```bash
--reads 'reads/*_R{1,2}.fastq.gz'      # always quote the glob
```

For Illumina-style names such as `*_R1_001.fastq.gz`, use `'reads/*_R{1,2}_001.fastq.gz'`.

### 7.2 Assemblies (alternative entry)

To analyse existing assemblies (e.g. from NCBI or a previous run), pass `--assemblies` instead of `--reads`:

* Steps 1–6 are skipped.
* The assembly gate and every downstream step run as usual.
* The sample ID is the file name without `.fasta`/`.fa`/`.fna` and without a trailing `_contigs`.

```bash
--assemblies 'genomes/*.fasta'
```

## 8. Usage

**Full run from reads (*A. baumannii*):**

```bash
nextflow run main.nf -profile conda \
  --reads 'reads/*_R{1,2}.fastq.gz' \
  --phix db/phix.fa \
  --kraken_db db/k2_standard_08gb \
  --reference db/ATCC17978.fasta \
  --card_json db/card.json \
  --outdir results
```

**Another species, e.g. *Klebsiella pneumoniae*:**

```bash
nextflow run main.nf -profile conda \
  --reads 'reads/*_R{1,2}.fastq.gz' \
  --kraken_db db/k2_standard_08gb \
  --target_species "Klebsiella pneumoniae" \
  --prokka_genus Klebsiella --prokka_species pneumoniae \
  --reference db/HS11286.fasta --card_json db/card.json
```

**Existing assemblies (pangenome and AMR only):**

```bash
nextflow run main.nf -profile conda \
  --assemblies 'genomes/*.fasta' --card_json db/card.json
```

**Resume after an interruption, or after changing a downstream parameter:**

```bash
nextflow run main.nf -profile conda <same options> -resume
```

**Laptop or small workstation:** add `--max_cpus 4 --max_memory 12.GB`.

## 9. Parameters

### Input / output

| Parameter | Default | Description |
|---|---|---|
| `--reads` | – | Glob of paired FASTQ files. Required unless `--assemblies` is given. |
| `--assemblies` | – | Glob of FASTA assemblies. Skips read processing. |
| `--outdir` | `results` | Output directory. |

### Read processing and species gate

| Parameter | Default | Description |
|---|---|---|
| `--phix` | – | PhiX FASTA. Leave it out to skip PhiX removal. |
| `--trimmer` | `trim_galore` | Read trimmer: `trim_galore` or `fastp`. |
| `--trim_quality` | `20` | Trim Galore: Phred cut-off for 3′ quality trimming. |
| `--trim_min_length` | `36` | Trim Galore: drop pairs with a read shorter than this. |
| `--trim_galore_fastqc` | `true` | Trim Galore: run FastQC on the trimmed reads. |
| `--trim_galore_args` | `""` | Extra Trim Galore arguments, passed through unchanged. |
| `--fastp_args` | `""` | Extra fastp arguments (only with `--trimmer fastp`). |
| `--kraken_db` | – | Kraken2/Bracken database folder. Required unless `--skip_species_check`. |
| `--read_len` | `100` | Bracken read length. Must match a `kmer_distrib` file in the database. |
| `--bracken_level` | `S` | Bracken taxonomic level. |
| `--target_species` | `Acinetobacter baumannii` | Expected top species. `""` checks purity only. |
| `--min_species_fraction` | `0.93` | Minimum Bracken fraction for the top species. |
| `--skip_species_check` | `false` | Skip Kraken2, Bracken and the gate. |

### Assembly

| Parameter | Default | Description |
|---|---|---|
| `--min_contig` | `500` | SKESA minimum contig length (bp). |
| `--max_contigs` | `200` | Assembly gate: maximum number of contigs. |

### Annotation, pangenome, comparative

| Parameter | Default | Description |
|---|---|---|
| `--prokka_genus` / `--prokka_species` | `Acinetobacter` / `baumannii` | Prokka organism. `--usegenus` is always set. |
| `--prokka_proteins` | – | Trusted protein FASTA for Prokka (optional). |
| `--prokka_centre` | `CDC` | Sequencing-centre tag used for `--compliant` IDs. |
| `--prokka_evalue` | `1e-6` | Prokka similarity e-value cut-off. |
| `--skip_pangenome` | `false` | Do not run Roary. |
| `--roary_identity` | `95` | Minimum BLASTP percentage identity. |
| `--roary_core` | `99` | Percentage of isolates a gene must be in to count as core. |
| `--roary_core_alignment` | `true` | Build the core-gene alignment (`-e`). |
| `--roary_mafft` | `false` | Use MAFFT instead of PRANK for the core alignment (faster). |
| `--reference` | – | Reference genome for FastANI, Snippy and QUAST. |
| `--card_json` | – | CARD `card.json` for RGI. |

### Resources

| Parameter | Default | Description |
|---|---|---|
| `--max_cpus` | the smaller of 16 and the available cores | Upper limit for any single task. |
| `--max_memory` | `16.GB` | Upper limit for any single task. |
| `--max_time` | `48.h` | Upper limit for any single task. |

## 10. Output

```
results/
├── 01_phix_removed/   <S>_phix.tsv                       (read counts; cleaned reads stay in work/)
├── 02_trim_galore/    *_trimming_report.txt, <S>_val_{1,2}_fastqc.{html,zip}
│   (or 02_fastp/      <S>_fastp.{html,json}  with --trimmer fastp)
├── 03_kraken2/        <S>.kraken2.report
├── 04_bracken/        <S>.bracken, <S>_bracken.report
├── 05_qc_gates/       species_check.tsv, assembly_qc.tsv
├── 06_assembly/       <S>_contigs.fasta
├── 07_fastani/        <S>.ani.txt
├── 08_prokka/<S>/     <S>.gff .gbk .faa .ffn .fna .tsv .txt …
├── 09_pangenome/roary/
│                      gene_presence_absence.{csv,Rtab}, summary_statistics.txt,
│                      core_gene_alignment.aln, accessory_binary_genes.fa.newick …
├── 10_snippy/         <S>/ (snps.vcf, snps.tab …)
│   └── core/          core.{aln,full.aln,tab,txt,vcf}
├── 11_rgi/            <S>.txt, <S>.json
├── 12_quast/quast/    report.{html,tsv} …
├── 13_multiqc/        multiqc_report.html, multiqc_data/
└── pipeline_info/     report.html, timeline.html, trace.txt, dag.html
```

## 11. Interpreting the results

**`05_qc_gates/*.tsv`.** Start here. Every sample missing from the downstream folders is listed with the reason.

**FastANI.**
* ANI of at least 95 % to the reference confirms the species.
* Values between 90 and 95 % within the *Acinetobacter* complex suggest *A. pittii*, *A. nosocomialis* or *A. seifertii* that culture misidentified.

**Roary.**
* `summary_statistics.txt` gives core, soft-core, shell and cloud gene counts. If the cloud fraction keeps growing as genomes are added, the pangenome is open.
* `gene_presence_absence.csv` can be loaded into Phandango, Scoary or R.

**Snippy `core/core.txt`.**
* It gives the aligned length and variant count for each sample. A sample with an unusually low `ALIGNED` share has probably diverged from the reference; consider a closer reference.
* `snp-dists core.aln` gives pairwise SNP distances. In the CRAB study [1], small SNP distances were used to define local genomic micro-clusters.

**RGI.**
* Focus on *Strict*/*Perfect* hits in the carbapenemase families: OXA-23/24/58-like, NDM and KPC.
* The intrinsic *A. baumannii* OXA-51-like and ADC (*ampC*) genes are expected in every isolate. Whether they give a carbapenem or cephalosporin phenotype depends on an upstream IS*Aba1* insertion, and RGI does **not** report that.

## 12. Computational requirements

Typical cost per sample, for a genome of about 4 Mb sequenced at about 100× with PE150 reads:

| Step | Typical cost |
|---|---|
| Kraken2 | Memory about equal to the database size (8 GB for Standard-8); a few minutes |
| SKESA | 4–8 GB, 5–15 min on 8 cores |
| Prokka | ~5–10 min on 8 cores |
| Roary + PRANK | The bottleneck: hours for >100 genomes. Use `--roary_mafft` to speed it up. |
| RGI | ~2–5 min |

The defaults fit a 16 GB workstation. For a large Kraken2 database (e.g. PlusPF, ~70 GB), either:

* raise `--max_memory` and the `KRAKEN2` memory in `nextflow.config`, or
* use the capped Standard-8/16 indexes.

## 13. Software versions

Versions pinned in `nextflow.config` (Bioconda):

| Tool | Version | Step |
|---|---|---|
| BWA | 0.7.18 | PhiX removal |
| SAMtools | 1.20 | PhiX removal |
| Trim Galore | 0.6.10 | read trimming (default) |
| fastp | 0.23.4 | read trimming (`--trimmer fastp`) |
| Kraken2 | 2.1.3 | classification |
| Bracken | 2.9 | abundance |
| SKESA | 2.5.1 | assembly |
| FastANI | 1.34 | ANI |
| Prokka | 1.14.6 | annotation |
| Roary | 3.13.0 | pangenome |
| Snippy | 4.6.0 | core SNPs |
| RGI | 6.0.3 | AMR |
| QUAST | 5.2.0 | assembly QC |
| MultiQC | 1.25 | report |

RGI results depend on **both** the RGI version and the CARD data release. Report the CARD version in any publication; `rgi database --version --local` shows it.

## 14. Testing

`test_data/` holds a tiny synthetic dataset: three read pairs, a PhiX stub, a reference and a CARD stub. It exercises every branch of the workflow without running the real tools:

```bash
nextflow run main.nf -profile test -stub
```

In a few seconds this checks:

* the channel logic
* both gates
* the guards on the pangenome and snippy-core steps
* output publishing

Real-tool validation on public *A. baumannii* data is planned. Until then, compare a few samples against earlier bash-based runs.

## 15. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `Process requirement exceeds available memory` | Lower `--max_memory` to what the machine has (e.g. `12.GB`). |
| Bracken: `kmer distribution file not found` | `--read_len` does not match a `database*mers.kmer_distrib` file in `--kraken_db`. |
| Every sample FAILS the species gate | Check the exact species name in a `.bracken` file and use it for `--target_species`. Also check `--min_species_fraction`. |
| `Pangenome skipped` warning | Fewer than 2 genomes passed both gates. |
| FastANI/Snippy output missing | `--reference` was not given. |
| RGI output missing | `--card_json` was not given. |
| Roary very slow | Add `--roary_mafft`, or use `--roary_core_alignment false` if only presence/absence is needed. |
| `--phix null` still runs PhiX removal | Leave the option out entirely instead. |
| Conda solving hangs | Check that mamba is installed, delete `work/conda/` and retry. |

## 16. Known limitations and roadmap

**v1.0.0 status:**

* The workflow logic has been verified with stub runs.
* Both QC gates passed real-mode tests on synthetic inputs.
* A full real-data run of the consolidated workflow has **not** been done yet.
* The Docker profile has not been tested.

Planned for v1.1. These are left out of v1.0.0 so that it stays faithful to the published analysis.

* MLST (Pasteur and Oxford schemes) and Kaptive K/OC-locus typing for *A. baumannii*
* CheckM2 completeness and contamination as an additional assembly gate
* Panaroo as an alternative pangenome tool that copes better with annotation errors
* IQ-TREE maximum-likelihood phylogeny and `snp-dists` on the core alignment
* Optional Gubbins recombination masking
* Detection of IS*Aba1* upstream of *bla*<sub>OXA-51-like</sub> and *bla*<sub>ADC</sub>
* Samplesheet (CSV) input

## 17. Citation

If you use this workflow, please cite the software (see `CITATION.cff`; the Zenodo DOI will be added on release) **and** the tools listed in [§18](#18-references).

> Algarni AT. *bacterial-wgs-nf: a Nextflow workflow for short-read bacterial genome assembly, pangenome and AMR analysis.* GitHub: https://github.com/atalgarni/bacterial-wgs-nf

## 18. References

1. Alquait AA, …, Algarni AT, *et al.* Phenotypic, PCR, and whole-genome sequencing characterization of carbapenem-resistant *Acinetobacter baumannii* in Makkah and Riyadh, Saudi Arabia: AST discordance and local genomic micro-clusters. *Front Microbiol.* 2026;17:1845444.
2. Li H, Durbin R. Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics.* 2009;25(14):1754–1760. doi:10.1093/bioinformatics/btp324
3. Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv.* 2013;1303.3997.
4. Danecek P, Bonfield JK, Liddle J, *et al.* Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. doi:10.1093/gigascience/giab008
5. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics.* 2018;34(17):i884–i890. doi:10.1093/bioinformatics/bty560
6. Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biol.* 2019;20:257. doi:10.1186/s13059-019-1891-0
7. Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. *PeerJ Comput Sci.* 2017;3:e104. doi:10.7717/peerj-cs.104
8. Souvorov A, Agarwala R, Lipman DJ. SKESA: strategic k-mer extension for scrupulous assemblies. *Genome Biol.* 2018;19:153. doi:10.1186/s13059-018-1540-z
9. Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nat Commun.* 2018;9:5114. doi:10.1038/s41467-018-07641-9
10. Seemann T. Prokka: rapid prokaryotic genome annotation. *Bioinformatics.* 2014;30(14):2068–2069. doi:10.1093/bioinformatics/btu153
11. Hyatt D, Chen GL, LoCascio PF, *et al.* Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics.* 2010;11:119. doi:10.1186/1471-2105-11-119
12. Laslett D, Canback B. ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. *Nucleic Acids Res.* 2004;32(1):11–16. doi:10.1093/nar/gkh152
13. Nawrocki EP, Eddy SR. Infernal 1.1: 100-fold faster RNA homology searches. *Bioinformatics.* 2013;29(22):2933–2935. doi:10.1093/bioinformatics/btt509
14. Kalvari I, Nawrocki EP, Ontiveros-Palacios N, *et al.* Rfam 14: expanded coverage of metagenomic, viral and microRNA families. *Nucleic Acids Res.* 2021;49(D1):D192–D200. doi:10.1093/nar/gkaa1047
15. Camacho C, Coulouris G, Avagyan V, *et al.* BLAST+: architecture and applications. *BMC Bioinformatics.* 2009;10:421. doi:10.1186/1471-2105-10-421
16. Page AJ, Cummins CA, Hunt M, *et al.* Roary: rapid large-scale prokaryote pan genome analysis. *Bioinformatics.* 2015;31(22):3691–3693. doi:10.1093/bioinformatics/btv421
17. Löytynoja A. Phylogeny-aware alignment with PRANK. *Methods Mol Biol.* 2014;1079:155–170. doi:10.1007/978-1-62703-646-7_10
18. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol Biol Evol.* 2013;30(4):772–780. doi:10.1093/molbev/mst010
19. Seemann T. Snippy: fast bacterial variant calling from NGS reads. 2015. https://github.com/tseemann/snippy
20. Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. *arXiv.* 2012;1207.3907.
21. Alcock BP, Huynh W, Chalil R, *et al.* CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res.* 2023;51(D1):D690–D699. doi:10.1093/nar/gkac920
22. Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. *Nat Methods.* 2015;12:59–60. doi:10.1038/nmeth.3176
23. Gurevich A, Saveliev V, Vyahhi N, Tesler G. QUAST: quality assessment tool for genome assemblies. *Bioinformatics.* 2013;29(8):1072–1075. doi:10.1093/bioinformatics/btt086
24. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics.* 2016;32(19):3047–3048. doi:10.1093/bioinformatics/btw354
25. Grüning B, Dale R, Sjödin A, *et al.* Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nat Methods.* 2018;15:475–476. doi:10.1038/s41592-018-0046-7
26. Di Tommaso P, Chatzou M, Floden EW, *et al.* Nextflow enables reproducible computational workflows. *Nat Biotechnol.* 2017;35:316–319. doi:10.1038/nbt.3820
27. Krueger F, James F, Ewels P, Afyounian E, Weinstein M, Schuster-Boeckler B, Hulselmans G, sclamons. FelixKrueger/TrimGalore. Zenodo. doi:10.5281/zenodo.5127898
28. Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal.* 2011;17(1):10–12. doi:10.14806/ej.17.1.200
29. Andrews S. FastQC: a quality control tool for high throughput sequence data. 2010. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

---

**Author:** Abdulmalek T. Algarni ([ORCID 0000-0002-2872-3724](https://orcid.org/0000-0002-2872-3724)). Released under the MIT License.
