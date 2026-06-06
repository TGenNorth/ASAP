.. |copy|   unicode:: U+000A9 .. COPYRIGHT SIGN

.. code-block:: none

    █████╗  ███████╗  █████╗  ██████╗
   ██╔══██╗ ██╔════╝ ██╔══██╗ ██╔══██╗
   ███████║ ███████╗ ███████║ ██████╔╝       
   ██╔══██║ ╚════██║ ██╔══██║ ██╔═══╝       
   ██║  ██║ ███████║ ██║  ██║ ██║
   ╚═╝  ╚═╝ ╚══════╝ ╚═╝  ╚═╝ ╚═╝
   ──────────────────────────────────────────────────────────────────────────────────────────
   Amplicon Sequencing Analysis Pipeline (ASAP)
   ──────────────────────────────────────────────────────────────────────────────────────────

**Pre-release** | Authors: Darrin Lemmer, W. Tanner Porter, *et al.*,
Pathogen & Microbiome Division, Translational Genomics Research Institute (TGen) |
License: |copy| TGen North (non-commercial)

**ASAP** is a start-to-finish Nextflow pipeline for targeted amplicon sequencing analysis.
The pipeline ingests demultiplexed sequencing reads from multiple platforms (Illumina,
Oxford Nanopore, PacBio), performs read-level quality control, reference-based alignment,
and produces quality-control metrics, SNPs, iSNVs, and consensus sequences. Outputs are
viewable through built-in reporting workflows or accessible as structured data for
project-specific downstream analysis. ASAP provides a flexible, open-source platform
suitable for users with limited bioinformatics experience, while offering advanced
customization for expert users. Previously, ASAP has been applied to pathogen and
antimicrobial resistance (AMR) detection and viral whole-genome assembly.

**Intended Uses:**

1. Pathogen detection and AMR profiling from targeted amplicon panels
2. Viral whole-genome assembly from tiled amplicon sequencing (e.g., SARS-CoV-2, RSV)
3. Multi-reference analysis for species differentiation or cross-panel quality control
4. High-resolution iSNV detection in mixed infections or heteroresistant populations

----

Overview
========

Targeted sequencing has become an essential tool across clinical research, ecology, and
infectious disease surveillance. By amplifying targeted genomic regions, this approach
provides high coverage depth and increased sensitivity in complex samples at a fraction of
the cost of non-targeted approaches. Target-specific primers and probes can be designed to
amplify highly specific genomic regions, enabling high taxonomic resolution, resistance
marker detection, virulence factor identification, or tiled coverage across a gene or
full genome.

ASAP is designed for **target-specific amplicon analysis** — aligning reads against one or
more references and using the resulting alignment files to synthesize final datasets
including amplicon presence/absence, coverage statistics, consensus sequences,
single-nucleotide polymorphism (SNP) tables, and intra-host single-nucleotide variants
(iSNVs).

**Highlights:**

- Supports Illumina short reads (paired-end and single-end), Oxford Nanopore, and PacBio
- Multi-reference analysis: multiple amplicons or species analyzed in a single run
- Primer masking, percent-identity filtering, and SMOR error correction for high-resolution variant calling
- Optional iVAR integration for primer trimming, variant calling, and consensus generation
- R-based post-processing (ASAP Tools) for SNP tables, coverage tables, QC figures, and FASTA export
- Full QC reporting via FastQC and MultiQC
- HPC-ready via SLURM with exponential retry and configurable resource scaling

----

Pipeline Summary
================

.. code-block:: none

   Reference Input (FASTA / GenBank / Excel / JSON)
             │
             ▼
   PREPARE_ASAP_JSON ──► GENERATE_REFERENCE_FASTA ──────────────┐
                                                                │
   Read Files (FASTQ) ──────────────────────────────────────────┤
             │                                                  │
             ▼                                                  │
   FastQC (initial QC)                                          │
             │                                                  │
             ▼                                                  │
   ┌─────────────────────────┐                                  │
   │ Illumina:  fastp        │                                  │
   │ ONT/PacBio: fastplong   │                                  │
   └─────────────────────────┘                                  │
             │                                                  │
   FastQC (post-trim QC)                                        │
             │                                                  │
             ▼                                                  ▼
   ┌──────────────────────────────────────────────────────────────┐
   │                      Alignment                               │
   │   Illumina:     Bowtie2 (default) or BWA-MEM                 │
   │   ONT/PacBio:   minimap2 (auto-selected)                     │
   └──────────────────────────────────────────────────────────────┘
             │
             ▼
   [Optional] Primer Masking        ← maskPrimers.py + BED file
             │
   [Optional] Percent-Identity Filtering  ← identityFilter.py
             │
   [Optional] SMOR Masking / SMOR Correction ← generateSMORbam.py
             │
             ├─────────────────────────────────────────────────────┐
             ▼                                                     ▼
   ┌──────────────────────────┐                    ┌───────────────────────────┐
   │  ASAP BAM Processing     │                    │  iVAR [optional]          │
   │  (newBamProcessor.py)    │                    │  ─ Primer trimming        │
   │  → per-sample XML        │                    │  ─ Variant calling (.tsv) │
   └──────────────────────────┘                    │  ─ Consensus FASTA        │
             │                                     └───────────────────────────┘
             ▼
   [Optional] ASAP Tools R Post-Processing
     ├── Per-sample XML → Rdata  (parallel)
     ├── Combine all Rdata       (gather)
     ├── Coverage Table          (Excel)
     ├── QC Figures              (HTML + JPG)
     ├── Consensus FASTA export
     ├── SNP → Amino Acid translation  (if GenBank provided)
     └── SNP / iSNV Table        (CSV ± Excel)
             │
   [Optional] OUTPUT_COMBINER + FORMAT_OUTPUT
             │  (combined XML → HTML report)
             ▼
   MultiQC (aggregated QC report)

----

Requirements
============

- **Nextflow** ≥ 23.04 (tested on 25.04.6) — must be available in your active environment
- **nf-schema** plugin 2.5.1 — loaded automatically via ``nextflow.config`` on first run
- **nf-test** ≥ 0.9.0 — required only to run the test suite
- **Singularity / Apptainer** (for containerized alignment and QC tools)
- **Conda / Mamba** (environments are built automatically from module YMLs — no manual setup required)
- A reference file in FASTA, GenBank, Excel (.xlsx), or JSON format

**Execution profiles:**

+---------------+---------------------+---------------------+-----------------------------------+
| Profile       | Executor            | Containers          | Best For                          |
+===============+=====================+=====================+===================================+
| ``slurm``     | SLURM (child jobs)  | Singularity + Conda | Production HPC runs               |
+---------------+---------------------+---------------------+-----------------------------------+
| ``conda``     | Local (current node)| Conda               | Interactive ``srun`` or laptop    |
+---------------+---------------------+---------------------+-----------------------------------+

**Key tool versions:**

+---------------+----------+--------------------------------------------------+
| Tool          | Version  | Purpose                                          |
+===============+==========+==================================================+
| fastp         | 0.23.4   | Illumina read QC and adapter trimming            |
+---------------+----------+--------------------------------------------------+
| fastplong     | 0.4.1    | ONT / PacBio read QC                            |
+---------------+----------+--------------------------------------------------+
| Bowtie2       | 2.x      | Short-read alignment                             |
+---------------+----------+--------------------------------------------------+
| BWA-MEM       | 0.7.x    | Short-read alignment (alternative)               |
+---------------+----------+--------------------------------------------------+
| minimap2      | 2.x      | Long-read alignment                              |
+---------------+----------+--------------------------------------------------+
| SAMtools      | 1.x      | BAM manipulation and indexing                    |
+---------------+----------+--------------------------------------------------+
| iVAR          | 1.4.4    | Primer trimming, variant calling, consensus      |
+---------------+----------+--------------------------------------------------+
| FastQC        | 0.12.1   | Per-sample read QC                               |
+---------------+----------+--------------------------------------------------+
| MultiQC       | latest   | Aggregated QC reporting                          |
+---------------+----------+--------------------------------------------------+

----

Installation
============

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/TGenNorth/ASAP.git
   cd ASAP/nextflow

   # Activate the conda environment that contains Nextflow (and nf-test for testing)
   conda activate <your-nextflow-env>

   # Verify Nextflow is available
   nextflow -version

   # View full parameter help
   nextflow run main.nf --help

Nextflow and nf-test must be available in your active environment — install them via
conda or follow the `Nextflow installation guide <https://www.nextflow.io/docs/latest/install.html>`_.
The **nf-schema** plugin (``nf-schema@2.5.1``) is declared in ``nextflow.config`` and
downloaded automatically on first run. Singularity containers and all Conda environments
for pipeline steps are also resolved automatically — no further manual setup is required.

----

Quick Start
===========

**Illumina paired-end reads with a GenBank reference (SLURM):**

.. code-block:: bash

   cd /path/to/ASAP/nextflow

   nextflow run main.nf \
     --read_dir        ./reads/ \
     --reference_input "./refs/*.gb" \
     --outdir          ASAP_Results \
     --file_name       MyRun \
     -profile slurm

**ONT long reads with minimap2:**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./ont_reads/ \
     --reference_input reference.gb \
     --outdir          ONT_Results \
     --file_name       ONT_Run \
     --technology      ont \
     -profile slurm

**JSON reference passed directly (no conversion):**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./reads/ \
     --reference_input assay.json \
     --outdir          Results \
     --file_name       MyRun \
     -profile slurm

**Resume a failed or interrupted run:**

.. code-block:: bash

   nextflow run main.nf [params] -resume

----

Reference Input Formats
=======================

ASAP accepts four reference formats via ``--reference_input``:

+----------+-----------------------+----------------------------------------------------------+
| Format   | Extension(s)          | Notes                                                    |
+==========+=======================+==========================================================+
| JSON     | ``.json``             | Used directly; no conversion step                        |
+----------+-----------------------+----------------------------------------------------------+
| GenBank  | ``.gb``, ``.gbk``    | Single or multiple files: ``"./refs/*.gb"``              |
+----------+-----------------------+----------------------------------------------------------+
| FASTA    | ``.fasta``, ``.fa``   | Single file; one amplicon entry per sequence             |
+----------+-----------------------+----------------------------------------------------------+
| Excel    | ``.xlsx``, ``.xls``   | Single file using the ASAP assay template                |
+----------+-----------------------+----------------------------------------------------------+

All non-JSON formats are converted to an internal JSON assay description by
``prepareJSONInput_nextflow.py`` before processing. The JSON encodes target names,
reference sequences, SNP positions of interest, and significance rules.

----

Pipeline Steps
==============

Step 1 — Quality Control
-------------------------

**Illumina reads** are processed by ``fastp`` [6]_, which performs adapter trimming,
quality filtering, and per-sample HTML/JSON reports.

**ONT and PacBio reads** are processed by ``fastplong`` [CITATION]_, a long-read
variant of fastp.

``FastQC`` [CITATION]_ runs on reads before and after trimming for per-sample QC
assessment.

+----------------------+-------------+----------------------------------------------------------+
| Parameter            | Default     | Description                                              |
+======================+=============+==========================================================+
| ``--technology``     | ``illumina``| Platform: ``illumina``, ``ont``, ``pacbio``              |
+----------------------+-------------+----------------------------------------------------------+
| ``--adapter_fasta``  | bundled     | Adapter FASTA; bundled Illumina adapters used by default |
+----------------------+-------------+----------------------------------------------------------+
| ``--fastp_extra_args``| ``""``     | Additional fastp flags (e.g. ``-l 100`` for min length)  |
+----------------------+-------------+----------------------------------------------------------+

Step 2 — Alignment
-------------------

Trimmed reads are aligned to ``reference.fasta`` (extracted from the assay JSON):

- **Illumina:** ``bowtie2`` (default) or ``bwa mem``
- **ONT / PacBio:** ``minimap2`` (selected automatically by ``--technology``)

All aligners produce a coordinate-sorted, indexed BAM published to
``sample_info/<sample>/bwa/``, ``bowtie2/``, or ``minimap2/`` respectively.
BWA and Bowtie2 also emit ``flagstat`` files for MultiQC.

+------------------------+------------+-------------------------------------------+
| Parameter              | Default    | Description                               |
+========================+============+===========================================+
| ``--aligner``          | ``bowtie2``| Aligner: ``bowtie2``, ``bwa``,``minimap2``|
+------------------------+------------+-------------------------------------------+
| ``--aligner_extra_args``| ``""``    | Additional arguments passed to the aligner|
+------------------------+------------+-------------------------------------------+

Step 3 — Primer Masking *(optional)*
--------------------------------------

Primer-derived base calls are masked in aligned reads to prevent them from inflating
or distorting SNP frequencies or consensus sequences. A BED-format TSV file specifying
primer coordinates is required.

For each read, if the alignment start (R1) or end (R2) falls within ``--wiggle`` bases
of a primer boundary, the primer region is masked: base quality scores are set to 0
and (by default) bases are replaced with ``N``. A per-amplicon log file tallies reads
with and without detected primer sequences, confirming correct masking.

+---------------------+---------+-------------------------------------------------------+
| Parameter           | Default | Description                                           |
+=====================+=========+=======================================================+
| ``--primer_file``   | ``null``| Path to primer BED file (required to enable masking)  |
+---------------------+---------+-------------------------------------------------------+
| ``--mask_primers``  | ``null``| Enable primer masking (auto-enabled when              |
|                     |         | ``--primer_file`` is provided; set ``false`` to force |
|                     |         | disable)                                              |
+---------------------+---------+-------------------------------------------------------+
| ``--wiggle``        | ``9``   | Bases outside primer boundary to include in mask      |
+---------------------+---------+-------------------------------------------------------+
| ``--mask_bam``      | ``true``| Replace masked bases with ``N`` in BAM sequence field |
+---------------------+---------+-------------------------------------------------------+
| ``--primer_only``   | ``false``| Retain only primer-overlapping reads; discard all    |
|                     |         | others after masking                                  |
+---------------------+---------+-------------------------------------------------------+

Step 4 — Percent-Identity Filtering *(optional)*
-------------------------------------------------

Read-level filtering removes reads that do not meet a minimum percent identity to their
aligned reference amplicon. Identity is computed via local Smith-Waterman alignment.

This is particularly valuable when near-neighbor organisms co-amplify with the target
and reads from the off-target organism must be excluded before variant calling (e.g.,
distinguishing *M. tuberculosis* from non-tuberculous mycobacteria).

+------------------+---------+---------------------------------------------------------------+
| Parameter        | Default | Description                                                   |
+==================+=========+===============================================================+
| ``--identity``   | ``null``| Minimum fractional identity threshold (e.g. ``0.97`` = 97%)  |
+------------------+---------+---------------------------------------------------------------+

Step 5 — SMOR Processing *(optional)*
--------------------------------------

Two complementary approaches leverage paired-end read overlap for Illumina error
correction, without requiring additional wet-lab steps.

**SMOR Masking (``--smor true``)**

Designed for assays where both reads in a pair are expected to fully overlap.
Non-overlapping read regions are masked, and discordant overlapping positions between
R1 and R2 are also masked. This approach yields an approximately 8-fold decrease in
overall error rate relative to unmasked data [CITATION]_, enabling resolution of minor
allele frequencies that would otherwise be indistinguishable from sequencing noise.

**SMOR Correction (``--smor_correction true``)**

Does not require full overlap. Within overlapping regions, discordant positions are
resolved by selecting the base with the higher Phred quality score (threshold: Q10
default). When reads agree, quality scores are combined to produce higher-confidence
calls. This approach is advantageous when reads partially overlap and quality degrades
toward the ends of R1 or R2.

+------------------------+---------+------------------------------------------+
| Parameter              | Default | Description                              |
+========================+=========+==========================================+
| ``--smor``             | ``false``| SMOR masking (full-overlap assays)      |
+------------------------+---------+------------------------------------------+
| ``--smor_correction``  | ``false``| SMOR correction (partial-overlap assays)|
+------------------------+---------+------------------------------------------+

Step 6 — ASAP BAM Processing
------------------------------

The core analysis step. ``newBamProcessor.py`` reads the assay JSON and the aligned
(optionally masked/filtered/SMOR'd) BAM to produce a per-sample XML containing:

- Aligned read counts per amplicon
- Per-position depth, breadth of coverage, and consensus sequence
- SNPs / iSNVs detected above the proportion and depth thresholds, with per-base distributions
- Regions of interest (ROI) sequences with DNA and amino acid translations
- Significance calls based on rules defined in the assay JSON

+--------------------------+----------+----------------------------------------------------------+
| Parameter                | Default  | Description                                              |
+==========================+==========+==========================================================+
| ``--depth``              | ``100``  | Minimum read depth to consider a position covered        |
+--------------------------+----------+----------------------------------------------------------+
| ``--breadth``            | ``0.8``  | Minimum breadth of coverage to call an amplicon present  |
+--------------------------+----------+----------------------------------------------------------+
| ``--proportion``         | ``0.1``  | Minimum allele frequency to call a SNP / iSNV            |
+--------------------------+----------+----------------------------------------------------------+
| ``--mutation_depth``     | ``5``    | Minimum read count to call a SNP / iSNV                  |
+--------------------------+----------+----------------------------------------------------------+
| ``--min_base_qual``      | ``5``    | Minimum Phred base quality score                         |
+--------------------------+----------+----------------------------------------------------------+
| ``--consensus_proportion``| ``0.8`` | Minimum frequency to call a consensus base (else ``N``)  |
+--------------------------+----------+----------------------------------------------------------+
| ``--fill_character``     | ``N``    | Character written at masked / gap positions (used by     |
|                          |          | SMOR masking and bam_processor)                          |
+--------------------------+----------+----------------------------------------------------------+
| ``--fill_gaps``          | ``n``    | Character written at zero-coverage positions in          |
|                          |          | consensus sequence                                       |
+--------------------------+----------+----------------------------------------------------------+
| ``--mark_deletions``     | ``_``    | Character written at deletion positions in consensus     |
+--------------------------+----------+----------------------------------------------------------+
| ``--whole_genome``       | ``false``| Skip per-sample consensus/depth arrays (WGS references)  |
+--------------------------+----------+----------------------------------------------------------+
| ``--asap_snps``          | ``true`` | Enable ASAP BAM processing (set ``false`` to skip)       |
+--------------------------+----------+----------------------------------------------------------+
| ``--combine_output``     | ``true`` | Combine per-sample XMLs and generate HTML report         |
+--------------------------+----------+----------------------------------------------------------+
| ``--stylesheet``         | bundled  | XSLT stylesheet for HTML report generation               |
+--------------------------+----------+----------------------------------------------------------+

Step 7 — iVAR Processing *(optional)*
---------------------------------------

iVAR [1]_ provides an alternative or complementary variant calling and consensus
generation workflow. Enable the full iVAR workflow with ``--ivar true``, or enable
individual steps independently.

+-------------------------------+---------+------------------------------------------------+
| Parameter                     | Default | Description                                    |
+===============================+=========+================================================+
| ``--ivar``                    | ``false``| Enable all iVAR steps (trim + variants + consensus) |
+-------------------------------+---------+------------------------------------------------+
| ``--ivar_trim``               | ``false``| iVAR primer trimming only (requires primer BED)|
+-------------------------------+---------+------------------------------------------------+
| ``--ivar_variants``           | ``false``| iVAR variant calling only                     |
+-------------------------------+---------+------------------------------------------------+
| ``--ivar_consensus``          | ``false``| iVAR consensus calling only                   |
+-------------------------------+---------+------------------------------------------------+
| ``--ivar_trim_extra_args``    | ``""``  | Additional ``ivar trim`` arguments              |
+-------------------------------+---------+------------------------------------------------+
| ``--ivar_variants_extra_args``| ``""``  | Additional ``ivar variants`` arguments          |
+-------------------------------+---------+------------------------------------------------+
| ``--ivar_consensus_extra_args``| ``""`` | Additional ``ivar consensus`` arguments         |
+-------------------------------+---------+------------------------------------------------+

Step 8 — ASAP Tools R Post-Processing *(optional)*
----------------------------------------------------

A suite of R scripts transforms per-sample XML outputs into tabular summaries,
figures, and FASTA files. Processing follows a fan-out / gather pattern:

.. code-block:: none

   PROCESS_XML_R         (per sample, parallel)  →  sample_info/<id>/rdata/
         │
         ▼
   PROCESS_COMBINE_RDATA (gather all)             →  sample_reports/rdata/ (Rdata)
         │                                           sample_reports/general_reports/ (CSV)
         ├── PROCESS_GENERATE_COV_TABLE           →  sample_reports/general_reports/ (Excel)
         ├── PROCESS_QC_PLOTS                     →  sample_reports/plots/ (HTML + JPG)
         ├── PROCESS_GENERATE_FASTA               →  sample_reports/fasta/ (FASTA)
         ├── PROCESS_SNPS_TO_AMINOACIDS           →  sample_reports/rdata/ (Rdata)
         │                                           sample_reports/snp_reports/ (CSV)
         └── PROCESS_GENERATE_SNP_TABLE           →  sample_reports/snp_reports/ (CSV ± Excel)

+--------------------------------------+----------+--------------------------------------------------+
| Parameter                            | Default  | Description                                      |
+======================================+==========+==================================================+
| ``--asaptools_processing``           | ``true`` | Enable R post-processing                         |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_cov_table``            | ``true`` | Generate coverage depth table                    |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_qc_plots``             | ``true`` | Generate QC figures                              |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_generate_fasta``       | ``true`` | Export consensus FASTA files                     |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_snp_table``            | ``true`` | Generate SNP / iSNV table                        |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_snp_table_xls``        | ``false``| Also export SNP table as Excel                   |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_positions_of_interest``| ``null`` | CSV of genomic positions to annotate in outputs  |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_genbank_location``     | ``null`` | GenBank file for amino acid annotation           |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_min_location_depth``   | ``99``   | Minimum depth to report a position in tables     |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_snp_proportion``       | ``null`` | Override allele frequency threshold for SNP table|
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_max_sample_snp_count`` | ``50``   | Max SNPs per sample before flagging as noisy     |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_samples_to_remove``    | ``null`` | Sample IDs to exclude from combined outputs      |
+--------------------------------------+----------+--------------------------------------------------+
| ``--asaptools_breadth_threshold``    | ``null`` | Minimum breadth for FASTA export (uses ``--breadth`` if unset) |
+--------------------------------------+----------+--------------------------------------------------+

Step 9 — MultiQC
-----------------

MultiQC [CITATION]_ aggregates reports from FastQC (initial and post-trim),
fastp/fastplong JSON, and SAMtools flagstat into a single interactive HTML report,
published to ``<outdir>/sample_reports/multiqc/``.

----

Unique Functionality
====================

Multi-Reference Analysis
-------------------------

ASAP is designed to analyze multiple amplicon targets simultaneously in a single run.
This supports:

- **Species co-analysis** — e.g., RSV-A and RSV-B from the same sequencing run
- **Multi-gene AMR panels** — multiple resistance genes across *M. tuberculosis*
- **Cross-contamination QC** — include neighboring-species references to detect and quantify cross-contamination between samples

Multi-reference analysis can be enabled by:

- Providing multiple sequences in a multi-FASTA or JSON reference file
- Specifying multiple GenBank files (``"./refs/*.gb"``)
- Using the ASAP Excel template with multiple amplicon rows
- Subsetting output annotation to specific positions via ``--asaptools_positions_of_interest``

Primer Masking
---------------

Primer masking prevents primer-derived base calls from distorting SNP frequencies or
consensus sequences — critical for amplicon assays where primers are included in
the sequenced region. A per-amplicon log file tallies reads with and without detected
primer sequences, enabling confirmation that the correct primer regions are being masked.
See `Step 3`_ for full parameter details.

Percent-Identity Filtering
---------------------------

Percent-identity filtering is especially valuable for:

- **Near-neighbor co-infections** — e.g., distinguishing *M. tuberculosis* from
  non-tuberculous mycobacteria when primers amplify both
- **Resistance gene specificity** — ensuring only reads from the precise target gene
  are used for resistance calling, preventing false calls from paralogs or related genes

See `Step 4`_ for full parameter details.

SMOR Masking and SMOR Correction
----------------------------------

SMOR approaches leverage paired-end read overlap to reduce sequencing errors without
requiring additional wet-lab steps. The ~8-fold error reduction achieved by SMOR Masking
enables resolution of minor allele frequencies (iSNVs) that would otherwise be
indistinguishable from sequencing noise, making it particularly powerful for:

- High-resolution AMR profiling in heteroresistant *M. tuberculosis* infections
- Detection of low-frequency viral variants in mixed infections
- Distinguishing true iSNVs from sequencing artifacts at low allele frequencies

See `Step 5`_ for the distinction between SMOR Masking and SMOR Correction.

ASAP Tools: Visualizations and Tabular Outputs
-----------------------------------------------

The R post-processing suite provides ready-to-use outputs for sample QC and downstream
analysis without requiring users to write custom analysis code:

**Coverage Table** — per-sample, per-amplicon statistics including breadth of coverage,
average depth, aligned reads, and depth at positions of interest. Depth thresholds flag
positions without adequate coverage, distinguishing true absence of SNPs from
insufficient data.

**QC Figures** — interactive HTML plots and static JPG figures showing depth of coverage
across samples, percent masked bases ("N"s), and SNP locations across the panel.

**SNP / iSNV Table** — a comprehensive table of detected variants across all samples and
amplicons. User-defined frequency and depth thresholds separate true iSNVs from noise.
When a GenBank reference is provided, coding-region SNPs are annotated with the
resulting amino acid change, enabling immediate identification of resistance-conferring
mutations.

**Consensus FASTA Export** — consensus sequences for each sample and amplicon, exported
at a user-defined breadth-of-coverage threshold. Suitable for downstream phylogenetic
analysis or genome assembly.

----

Output Structure
================

.. code-block:: none

   <outdir>/
   ├── pipeline_info/
   │   ├── sample_read_type_summary.tsv        # Sample IDs and read type (PE / SE)
   │   ├── dag_<name>_<timestamp>.png          # Pipeline DAG
   │   ├── report_<name>_<timestamp>.html      # Per-process execution report
   │   ├── trace_<name>_<timestamp>.txt        # Per-task resource usage
   │   └── timeline_<name>_<timestamp>.html    # Interactive job timeline
   │
   ├── reference/
   │   ├── reference.fasta                     # Extracted amplicon reference sequences
   │   └── bwa_index/ or bt2_index/            # Aligner index files
   │
   ├── sample_info/<sample>/
   │   ├── fastqc_initial/                     # Pre-trim FastQC reports
   │   ├── fastp/ or fastp_long/               # Trimmed reads + QC HTML / JSON
   │   ├── fastqc_post/                        # Post-trim FastQC reports
   │   ├── bwa/ or bowtie2/ or minimap2/       # Sorted, indexed BAM + flagstat
   │   ├── xml/
   │   │   └── <sample>.xml                    # Per-sample ASAP results
   │   ├── rdata/
   │   │   ├── <sample>_XML_Data.Rdata         # Per-sample R data object
   │   │   └── <sample>_Summary.csv            # Per-sample summary table
   │   ├── mask_primers/                       # Primer-masked BAM + masking log/TSV
   │   ├── identity_filter/                    # Identity-filtered BAM + log
   │   ├── smor/ or smor_correction/           # SMOR-processed BAM + log
   │   ├── ivar_trim/                          # iVAR-trimmed BAM
   │   ├── ivar_variants/                      # iVAR variant TSV
   │   └── ivar_consensus/                     # iVAR consensus FASTA
   │
   └── sample_reports/
       ├── <name>_analysis.xml                 # Combined XML (all samples)
       ├── <name>_report.html                  # Combined HTML report
       ├── multiqc/
       │   └── <name>_multiqc.html             # Aggregated QC report
       ├── rdata/
       │   ├── <name>_ASAP_Data.Rdata          # Combined R data (all samples)
       │   └── SNP_Amino_Acid_Table.Rdata      # Amino acid changes (GenBank required)
       ├── general_reports/
       │   ├── <name>_Summary.csv              # Combined summary table
       │   └── <name>_coverage_table.xlsx      # Coverage depth table
       ├── plots/
       │   ├── *.html                          # QC figures (interactive)
       │   └── *.jpg                           # QC figures (static)
       ├── snp_reports/
       │   ├── *.csv                           # SNP / iSNV table + amino acid changes
       │   └── *.xlsx                          # SNP table Excel (if --asaptools_snp_table_xls)
       └── fasta/
           └── *.fasta                         # Consensus FASTA exports

----

Example Commands
================

**TB amplicon panel — primer masking, identity filtering, SMOR correction, full asaptools:**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir          ./reads/ \
     --reference_input   ./refs/H37Rv.gb \
     --outdir            TB_Results \
     --file_name         TB_Run \
     --aligner           bwa \
     --primer_file       ./primers/TB_primers.bed \
     --mask_primers      true \
     --identity          0.97 \
     --smor_correction   true \
     --fastp_extra_args  "-l 100" \
     --proportion        0.01 \
     --asaptools_positions_of_interest ./genes/H37Rv_genes.csv \
     -profile slurm

**SARS-CoV-2 tiled amplicon — iVAR variant calling and consensus:**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./sc2_reads/ \
     --reference_input ./refs/SC2_Reference.json \
     --outdir          SC2_Results \
     --file_name       SC2_Run \
     --aligner         bwa \
     --primer_file     ./primers/SC2_primers.bed \
     --mask_primers    true \
     --ivar            true \
     --combine_output  true \
     -profile slurm

**RSV multi-GenBank reference (RSV-A + RSV-B simultaneously):**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./rsv_reads/ \
     --reference_input "./refs/RSV*.gb" \
     --outdir          RSV_Results \
     --file_name       RSV_Run \
     --aligner         bwa \
     --combine_output  true \
     -profile slurm

**ONT long reads with minimap2:**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./ont_reads/ \
     --reference_input ./refs/reference.gb \
     --outdir          ONT_Results \
     --file_name       ONT_Run \
     --technology      ont \
     --combine_output  true \
     -profile slurm

**ASAP Tools disabled — XML and HTML report only:**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir              ./reads/ \
     --reference_input       ./refs/reference.json \
     --outdir                Results \
     --file_name             Run \
     --asaptools_processing  false \
     --combine_output        true \
     -profile slurm

----

Performance Considerations
==========================

**SLURM resource scaling:** All processes use exponential retry. Memory and wall-time
double on each retry attempt (``memory = { N.GB * 2**(task.attempt-1) }``), preventing
transient resource limits from aborting runs.

**R post-processing memory:** ``PROCESS_COMBINE_RDATA`` (default 100 GB) and
``PROCESS_SNPS_TO_AMINOACIDS`` (default 30 GB) are the most memory-intensive steps,
as they load all per-sample Rdata objects simultaneously. For large runs (100+ samples),
ensure sufficient memory is available on the target SLURM partition.

**SMOR and identity filtering:** These steps add computational overhead but substantially
improve variant call quality. For routine screening, they can be omitted; for
high-resolution iSNV detection or resistance profiling, they are strongly recommended.

**Primer masking:** Runs on a single CPU and completes quickly relative to alignment.
The BAM re-sort step (``--mask_bam true``) adds a small overhead but ensures downstream
tools receive a correctly sorted BAM.

----

Running the Test Suite
======================

End-to-end tests are defined in ``tests/ASAP_EtE.nf.test`` and executed via
``nf-test``. A SLURM submission script is provided for convenience:

.. code-block:: bash

   cd /path/to/ASAP/nextflow
   mkdir -p logs/nf-test        # must exist before sbatch

   # Run all tests
   sbatch run_tests.sh

   # First run — generate snapshots for stable outputs
   sbatch run_tests.sh --update-snapshot

   # Run a single data type
   sbatch --job-name=ASAP_tb run_tests.sh --tag tb --keep-going

Available test tags: ``help``, ``rsv``, ``tb``, ``sc2``, ``bwa``, ``bowtie2``,
``minimap2``, ``multi_gb``, ``json_input``, ``excel_input``, ``fasta_input``,
``paired_end``, ``single_end``, ``ont``, ``primer_masking``, ``identity_filter``,
``smor``, ``ivar``, ``asaptools``, ``combine_output``, ``whole_genome``.

----

License
=======

Copyright |copy| The Translational Genomics Research Institute (TGen).
See the included ``LICENSE`` document.

Available for academic and research use under a license from TGen that is free for
non-commercial use. Distributed on an "AS IS" basis without warranties or conditions
of any kind, either express or implied.

----

Contact
=======

| TGen North
| 3051 W Shamrell Blvd Ste 106
| Flagstaff, AZ 86001-9435

| Darrin Lemmer — dlemmer@tgen.org
| W. Tanner Porter — tporter@tgen.org

Issues and feature requests:
https://github.com/TGenNorth/ASAP/issues

----

References
==========

.. [1] Grubaugh ND, Gangavarapu K, Quick J, *et al.* An amplicon-based sequencing
       framework for accurately measuring intrahost virus diversity using PrimalSeq
       and iVar. *Genome Biology*. 2019;20(1):8.

.. [2] Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2.
       *Nature Methods*. 2012;9:357–359.

.. [3] Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler
       Aligner. *Bioinformatics*. 2009;25(14):1754–1760.

.. [4] Li H. Minimap2: pairwise alignment for nucleotide sequences.
       *Bioinformatics*. 2018;34(18):3094–3100.

.. [5] Di Tommaso P, Chatzou M, Floden EW, *et al.* Nextflow enables reproducible
       computational workflows. *Nature Biotechnology*. 2017;35(4):316–319.

.. [6] Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ
       preprocessor. *Bioinformatics*. 2018;34(17):i884–i890.

.. [7] Li H, Handsaker B, Wysoker A, *et al.* The Sequence Alignment/Map Format
       and SAMtools. *Bioinformatics*. 2009;25(16):2078–2079.

.. [8] Andrews S. FastQC: A quality control tool for high throughput sequence data.
       2010. http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

.. [9] Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results
       for multiple tools and samples in a single report.
       *Bioinformatics*. 2016;32(19):3047–3048.

.. [10] [SMOR citation — TGen internal or forthcoming publication]

.. [CITATION] fastplong — citation pending.

.. [CITATION] FastQC — Andrews S. 2010 (see [8]_).

.. [CITATION] MultiQC — Ewels P. *et al.* 2016 (see [9]_).

----

AI Development Assistance
==========================

Portions of this pipeline — including module logic, subworkflow design, parameter
schema, test suite, and documentation — were developed with assistance from
**Claude Sonnet 4.6** (Anthropic). AI assistance was used as a coding and design
collaborator under active human direction and review. All scientific decisions,
parameter choices, and pipeline architecture reflect the work of the authors.
