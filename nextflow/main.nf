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
include { IVAR_TRIM           } from './modules/ivar/trim/'
include { IVAR_VARIANTS           } from './modules/ivar/variants/'
include { IVAR_CONSENSUS           } from './modules/ivar/consensus/'
include { MULTIQC                 } from './modules/multiqc'

// Call validation in the global scope.
validateParameters()

workflow {
    
    log.info paramsSummaryLog(workflow)

    // Generate reference fasta
    def assay_json = file(params.json).toAbsolutePath()
    json_ch = Channel.value(assay_json)
    def ref_fasta = GENERATE_REFERENCE_FASTA(json_ch)

    // Find and pair the read files
    def reads_dir_abs = file(params.read_dir).toAbsolutePath()
    def pattern = "${reads_dir_abs}/*_R{1,2}_001.fastq.gz"
    Channel
        .fromFilePairs(pattern)
        .ifEmpty { error "No paired-end FASTQ files matched in: ${reads_dir_abs}" }
        .set { paired_reads }
    
// --- STEP 1: FastQC Initial ---
    paired_reads
        .map { id, reads -> [ [id: id, status: 'initial'], reads ] } 
        .set { ch_for_fastqc_initial } 
    
    FASTQC_INITIAL(ch_for_fastqc_initial)

// --- STEP 2: Fastp Trimming ---
    def default_adapter_path = "${baseDir}/../asap/illumina_adapters_all.fasta"
    def adapter_path = params.adapter_fasta ?: default_adapter_path
    def adapter_file_obj = file(adapter_path, checkIfExists: true)

    adapter_fasta_ch = Channel.value(adapter_file_obj)

    def fastp_out = RUN_FASTP(paired_reads, adapter_fasta_ch)

    // Re-structure using .trimmed_reads to get the 5-item tuple
    fastp_out.trimmed_reads
        .map { id, r1, r2, html, json -> [ [id: id, status: 'post_process'], [r1, r2] ] }
        .set { ch_for_fastqc_post }

    // Re-structure for the aligner
    def trimmed_reads = fastp_out.trimmed_reads.map { id, r1, r2, html, json -> tuple(id, r1, r2) }

// --- STEP 3: Rerun Fastqc --- 
    FASTQC_POST(ch_for_fastqc_post)

// --- STEP 4: Align reads ---
    def aligned_bams
    def ch_flagstats = Channel.empty() 

    switch(params.aligner.toLowerCase()) {
        case 'bowtie2':
            def index_dir = BUILD_BOWTIE2_INDEX(ref_fasta)
            ALIGN_BOWTIE2(trimmed_reads.combine(index_dir))
            aligned_bams = ALIGN_BOWTIE2.out.bam_output
            ch_flagstats = ALIGN_BOWTIE2.out.flagstat 
            break
        case 'bwa':
            def index_dir = BUILD_BWA_INDEX(ref_fasta)
            ALIGN_BWA(trimmed_reads.combine(index_dir))
            aligned_bams = ALIGN_BWA.out.bam_output
            ch_flagstats = ALIGN_BWA.out.flagstat 
            break
        default:
            error "Unknown aligner: ${params.aligner}. Use 'bwa' or 'bowtie2'."
    }

// --- Optional STEP 5: Primer masking ---
    def primer_bed_path = params.primer_file ? file(params.primer_file).toAbsolutePath() : null
    def perform_masking = params.mask_primers || (params.primer_file && params.mask_primers != false)

    if(perform_masking) {
        if (!primer_bed_path) {
            error "ERROR: '--mask_primers' is enabled, but no BED file was provided via '--primer_file'."
        }
        MASK_PRIMERS(aligned_bams.combine(Channel.value(primer_bed_path)))
        aligned_bams = MASK_PRIMERS.out.mask_primers_output
    }

// --- Optional STEP 6: Identity Filtering ---
    if(params.identity) {
        (aligned_bams, identity_filter_logging) = IDENTITY_FILTER(aligned_bams)
    }

// --- Optional STEP 7: SMOR Processing ---
    if(params.smor) {
        (aligned_bams, smor_logging) = SMOR(aligned_bams)
    }

// --- STEP 8 & 10: Forking the BAM channel ---
    // Ensure keys match: it[0] is the String ID
    aligned_bams.multiMap { it ->
        asap: [ it[0], it[1], it[2] ]   
        ivar: [ [id: it[0]], it[1] ]    
    }.set { ch_split }
    
    // Crucial: Map the BAI channel to use the same [id: ID] map key as ch_split.ivar
    def bam_indices = aligned_bams.map{ it -> [ [id: it[0]], it[2] ] }
    
    if(params.asap_snps) {
        def xml_output = PROCESS_BAM(ch_split.asap.combine(json_ch))

        if(params.combine_output) {
            def xmls = xml_output.map { id, f -> f }.collect()
            def final_xml = OUTPUT_COMBINER(xmls)
            def default_stylesheet = file("${workflow.projectDir}/../output_transforms/ASAP_fulldetails_web.xsl")
            def stylesheet_path = params.stylesheet ? file(params.stylesheet) : default_stylesheet
            FORMAT_OUTPUT(final_xml, Channel.value(stylesheet_path))
        }
    }

// --- Optional STEP 10: iVAR Trimming ---
    ch_bam_for_ivar = Channel.empty()

    if (params.ivar || params.ivar_trim) {
        if (!params.primer_file) {
            error "ERROR: iVar trimming requested, but no primer file provided via '--primer_file'."
        }
        // join() pairs [id:X, bam] with [id:X, bai] -> [id:X, bam, bai]
        IVAR_TRIM (
            ch_split.ivar.join(bam_indices), 
            primer_bed_path
        )
        ch_bam_for_ivar = IVAR_TRIM.out.bam
    } else {
        ch_bam_for_ivar = ch_split.ivar
    }

// --- Optional STEP 11 & 12: iVAR Variants & Consensus --- 
    if (params.ivar || params.ivar_variants) {
        IVAR_VARIANTS (ch_bam_for_ivar, ref_fasta, true)
    }
    
    if (params.ivar || params.ivar_consensus) {
        IVAR_CONSENSUS (ch_bam_for_ivar, ref_fasta, true)
    }

// --- STEP 13: MultiQC Integration ---
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_INITIAL.out.zip.map{ it[1] }.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.map{ it[1] }.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(RUN_FASTP.out.json.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_flagstats.collect().ifEmpty([]))

    def ivar_v_stats = (params.ivar || params.ivar_variants) ? 
        IVAR_VARIANTS.out.tsv.map{ meta, tsv -> tsv }.collect().ifEmpty([]) : 
        Channel.empty()

    def ivar_t_stats = (params.ivar || params.ivar_trim) ? 
        IVAR_TRIM.out.log.map{ meta, log -> log }.collect().ifEmpty([]) : 
        Channel.empty()

    MULTIQC (
        ch_multiqc_files
            .mix(ivar_v_stats)
            .mix(ivar_t_stats)
            .collect(),
        [], [], [], [], []
    )
}

// Sub-workflows for organized FastQC calling
workflow FASTQC_INITIAL {
    take: reads
    main: FASTQC(reads)
    emit: zip = FASTQC.out.zip
}

workflow FASTQC_POST {
    take: reads
    main: FASTQC(reads)
    emit: zip = FASTQC.out.zip
}