#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// Import nf-schema functions
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// Import all processes
include { FASTQC } from './modules/fastqc'
include { RUN_FASTP } from './modules/fastp'
include { 
    BUILD_BWA_INDEX
    ALIGN_BWA
} from './modules/bwa'
include { 
    BUILD_BOWTIE2_INDEX
    ALIGN_BOWTIE2
} from './modules/bowtie2'
include {
    GENERATE_REFERENCE_FASTA
    MASK_PRIMERS
    IDENTITY_FILTER
    SMOR
    PROCESS_BAM
    OUTPUT_COMBINER
    FORMAT_OUTPUT
} from './modules/asap'
include { IVAR_VARIANTS           } from './modules/ivar/variants/'

// Call validation in the global scope. The plugin handles the --help flag.
validateParameters()

workflow {
    
    // Everything else goes inside the workflow block
    log.info paramsSummaryLog(workflow)

    // Grap the assay description json file and genereate reference fasta
    def assay_json = file(params.json).toAbsolutePath()
    json_ch = Channel.value(assay_json)
    def ref_fasta = GENERATE_REFERENCE_FASTA(json_ch)

    // Find and pair the read files for processing
    def reads_dir_abs = file(params.read_dir).toAbsolutePath()
    def pattern = "${reads_dir_abs}/*_R{1,2}_001.fastq.gz"
    Channel
        .fromFilePairs(pattern)
        .ifEmpty { error "No paired-end FASTQ files matched in: ${reads_dir_abs}" }
        .set { paired_reads }
    
// --- STEP 1: FastQC Initial ---
    // Run FastQC on reads.
    paired_reads
        .map { id, reads -> [ [id: id, status: 'initial'], reads ] } // Added status here
        .set { ch_for_fastqc_initial } // Must match the variable name below
    
    FASTQC_INITIAL(ch_for_fastqc_initial)

// --- STEP 2: Fastp Trimming ---
    // Get the adapter fasta to use for trimming
    def default_adapter = file("${workflow.projectDir}/../asap/illumina_adapters_all.fasta")
    def adapter_path = params.adapter_fasta ? file(params.adapter_fasta) : default_adapter
    Channel
        .value(adapter_path)
        .set { adapter_fasta }

    // Run fastp for adapter and quality trimming
    def fastp_out = RUN_FASTP(paired_reads, adapter_fasta)
    // Re-structure fastp output for the next FastQC run
    fastp_out
        .map { id, r1, r2, html, json -> [ [id: id, status: 'post_process'], [r1, r2] ] }
        .set { ch_for_fastqc_post }
    // Re-structure fastp output for the aligner
    def trimmed_reads = fastp_out.map { id, r1, r2, html, json -> tuple(id, r1, r2) }

// --- STEP 3: Rerun Fatqc ---

    FASTQC_POST(ch_for_fastqc_post)

// --- STEP 4: Align Fastp trimmed reads ---
    // Align reads using either bowtie2 or bwa
    def aligned_bams
    switch(params.aligner.toLowerCase()) {
        case 'bowtie2':
            def index_dir = BUILD_BOWTIE2_INDEX(ref_fasta)
            aligned_bams = ALIGN_BOWTIE2(trimmed_reads.combine(index_dir))
            break
        case 'bwa':
            def index_dir = BUILD_BWA_INDEX(ref_fasta)
            aligned_bams = ALIGN_BWA(trimmed_reads.combine(index_dir))
            break
        default:
            error "Unknown aligner: ${params.aligner}. Use 'bwa' or 'bowtie2'."
    }

// --- Optional STEP 5: Primer masking ---
    // Optionally run primer masking
    if(params.primer_file) {
        def mask_primers_logging
        Channel.of(file(params.primer_file).toAbsolutePath()).set { primer_file_ch }
        (aligned_bams,mask_primers_logging) = MASK_PRIMERS(aligned_bams.combine(primer_file_ch))
    }

// --- Optional STEP 6: Identity Filtering ---
    // Optionally run identity filtering
    if(params.identity) {
        def identity_filter_logging
        (aligned_bams, identity_filter_logging) = IDENTITY_FILTER(aligned_bams)
    }

// --- Optional STEP 7: SMOR Processing ---
    if(params.smor) {
        def smor_logging
        (aligned_bams, smor_logging) = SMOR(aligned_bams)
    }

// --- STEP 8 & 10: Forking the BAM channel ---
    // We split aligned_bams into two independent channels so they can both 
    // run their respective 6 samples without "stealing" from each other.
    aligned_bams.multiMap { it ->
        // it[0]=id, it[1]=bam, it[2]=bai
        asap: [ it[0], it[1], it[2] ]   // Pass ID, BAM, and BAI (3 items)
        ivar: [ [id: it[0]], it[1] ]    // Pass Meta and BAM for iVar
    }.set { ch_split }

// --- STEP 8: Bam Processor (SNP calling and XML creation) ---
    // We take the forked asap channel [id, bam, bai] and add json
    def xml_output = PROCESS_BAM(ch_split.asap.combine(json_ch))

// --- Optional STEP 9: Output Combiner --- 
    // Optionally run output combiner and transformation
    if(params.combine_output) {
        def xmls = xml_output.map { id, f -> f }.collect()
        def final_xml = OUTPUT_COMBINER(xmls)

        // Get the stylesheet to use for transformation
        def default_stylesheet = file("${workflow.projectDir}/../output_transforms/ASAP_fulldetails_web.xsl")
        def stylesheet_path = params.stylesheet ? file(params.stylesheet) : default_stylesheet
        Channel
            .value(stylesheet_path)
            .set { stylesheet_ch }
        def transformation = FORMAT_OUTPUT(final_xml, stylesheet_ch)
    }

// --- Optional STEP 10: iVAR Variant Calling --- 
    IVAR_VARIANTS (
        ch_split.ivar,
        ref_fasta,
        params.primer_file ? file(params.primer_file) : [],
        [], 
        true 
    )
}

workflow FASTQC_INITIAL {
    take: reads
    main: FASTQC(reads)
}

workflow FASTQC_POST {
    take: reads
    main: FASTQC(reads)
}