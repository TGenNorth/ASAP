#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// Import processes
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

workflow {

    // Grap the assay description json file and genereate reference fasta
    def assay_json = file(params.json).toAbsolutePath()
    Channel.of(assay_json).set { json_ch }
    def ref_fasta = GENERATE_REFERENCE_FASTA(json_ch)

    // Find and pair the read files for processing
    if (!params.read_dir) {
        error "Please specify the readfile location with: --read_dir <path>"
    }
    def reads_dir_abs = file(params.read_dir).toAbsolutePath()
    def pattern = "${reads_dir_abs}/*_R{1,2}_001.fastq.gz"
    Channel
        .fromFilePairs(pattern)
        .ifEmpty { error "No paired-end FASTQ files matched in: ${reads_dir_abs}" }
        .set { paired_reads }

    // Get the adapter fasta to use for trimming
    def default_adapter = file("${workflow.projectDir}/../asap/illumina_adapters_all.fasta")
    def adapter_path = params.adapter_fasta ? file(params.adapter_fasta) : default_adapter
    Channel
        .value(adapter_path)
        .set { adapter_fasta }

    // Run fastp for adapter and quality trimming
    def fastp_out = RUN_FASTP(paired_reads, adapter_fasta)
    def trimmed_reads = fastp_out.map { sample_id, r1, r2, html, json -> tuple(sample_id, r1, r2) }

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

    // Optionally run primer masking
    if(params.primer_file) {
        def mask_primers_logging
        Channel.of(file(params.primer_file).toAbsolutePath()).set { primer_file_ch }
        (aligned_bams,mask_primers_logging) = MASK_PRIMERS(aligned_bams.combine(primer_file_ch))
    }

    // Optionally run identity filtering
    if(params.identity) {
        def identity_filter_logging
        (aligned_bams, identity_filter_logging) = IDENTITY_FILTER(aligned_bams)
    }

    // Optionally run SMOR
    if(params.smor) {
        def smor_logging
        (aligned_bams, smor_logging) = SMOR(aligned_bams)
    }

    // Run bam processor
    def xml_output = PROCESS_BAM(aligned_bams.combine(json_ch))
    
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
}

