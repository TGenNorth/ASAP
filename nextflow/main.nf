#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// Import all processes
include { FASTQC } from './modules/fastqc'
include { RUN_FASTP } from './modules/fastp'
include { FASTPLONG } from './modules/fastplong'
include { BUILD_BWA_INDEX; ALIGN_BWA } from './modules/bwa'
include { BUILD_BOWTIE2_INDEX; ALIGN_BOWTIE2 } from './modules/bowtie2'
include { MINIMAP2_INDEX} from './modules/minimap2/index'
include { MINIMAP2_ALIGN } from './modules/minimap2/align'
include {
    GENERATE_REFERENCE_FASTA; MASK_PRIMERS; IDENTITY_FILTER; SMOR
    PROCESS_BAM; OUTPUT_COMBINER; FORMAT_OUTPUT
} from './modules/asap'
include { PROCESS_XML_R } from './modules/asap_tools'
include { IVAR_TRIM } from './modules/ivar/trim/'
include { IVAR_VARIANTS } from './modules/ivar/variants/'
include { IVAR_CONSENSUS } from './modules/ivar/consensus/'
include { MULTIQC } from './modules/multiqc'

validateParameters()

workflow {

    log.info paramsSummaryLog(workflow)
    
    // --- SETUP: Reference Generation ---
    def assay_json = file(params.json).toAbsolutePath()
    json_ch = Channel.value(assay_json)
    def ref_fasta = GENERATE_REFERENCE_FASTA(json_ch).ref_fasta
    
    // --- STEP 0: Input Handling ---
    def reads_dir_abs = file(params.read_dir).toAbsolutePath()
    def search_pattern = "${reads_dir_abs}/*.{fastq.gz,fq.gz}"

    def ch_raw_reads = Channel
        .fromFilePairs(search_pattern, checkIfExists: true, size: -1) { file -> 
            file.name.replaceAll(/(_R1|_R2|_[12])?(_001)?\.(fastq|fq)\.gz$/, '') 
        }
        .map { id, files ->
            def is_single = files.size() == 1
            [ [id: id, single_end: is_single], files ]
        }
    
    ch_raw_reads
        .map { meta, files -> 
            "${meta.id}\t${meta.single_end ? 'Single-End' : 'Paired-End'}\t${files.join(', ')}" 
        }
        .collectFile(
            name: 'sample_read_type_summary.tsv', 
            keepHeader: true, 
            newLine: true, 
            storeDir: "${params.outdir}/pipeline_info"
        ) { "Sample_ID\tType\tFiles" }

    def ch_raw_reads_for_pipeline = ch_raw_reads
    
    // --- STEP 1: FastQC Initial ---
    def ch_for_fastqc_initial = ch_raw_reads_for_pipeline
        .map { meta, reads -> [ meta.clone() << [status: 'initial'], reads ] } 
    
    FASTQC_INITIAL(ch_for_fastqc_initial)

    // --- STEP 2: Trimming Branch ---
    def ch_for_fastqc_post
    def ch_trimmed_for_align
    def ch_trim_json_for_multiqc

    if (params.technology.toLowerCase() == 'illumina') {
        def adapter_path = params.adapter_fasta ?: "${baseDir}/../asap/illumina_adapters_all.fasta"
        adapter_fasta_ch = Channel.value(file(adapter_path, checkIfExists: true))

        def fastp_out = RUN_FASTP(ch_raw_reads_for_pipeline, adapter_fasta_ch)

        ch_trimmed_for_align = fastp_out.trimmed_reads
        ch_trim_json_for_multiqc = fastp_out.json
        
        ch_for_fastqc_post = fastp_out.trimmed_reads
            .map { meta, reads -> [ meta.clone() << [status: 'post_process'], reads ] }

    } else if (params.technology.toLowerCase() == 'ont' || params.technology.toLowerCase() == 'ont.v14' || params.technology.toLowerCase() == 'pacbio') {
        
        def ch_longreads_unwrapped = ch_raw_reads_for_pipeline
            .map { meta, reads -> [ meta, reads[0] ] }

        def fastplong_out = FASTPLONG(ch_longreads_unwrapped)

        ch_trimmed_for_align = fastplong_out.reads
        ch_trim_json_for_multiqc = fastplong_out.json

        ch_for_fastqc_post = fastplong_out.reads
            .map { meta, reads -> [ meta.clone() << [status: 'post_process'], reads ] }

    } else {
        error "Unknown technology: ${params.technology}. Valid: illumina, ont, pacbio"
    }

    // --- STEP 3: Rerun Fastqc --- 
    FASTQC_POST(ch_for_fastqc_post)

    // --- STEP 4: Align reads ---
    def ch_aligned_with_meta
    def ch_flagstats = Channel.empty()

    def use_minimap2 = (
        params.technology.toLowerCase() == 'ont' || 
        params.technology.toLowerCase() == 'pacbio' || 
        params.aligner.toLowerCase() == 'minimap2'
    )

    if (use_minimap2) {
        def minimap_index = MINIMAP2_INDEX(ref_fasta.map{ it -> [[id: it.baseName], it] })
        MINIMAP2_ALIGN(
            ch_trimmed_for_align, 
            minimap_index.index.collect(), 
            true, "bai", false, true
        )
        ch_aligned_with_meta = MINIMAP2_ALIGN.out.bam_output
    } else {
        switch(params.aligner.toLowerCase()) {
            case 'bowtie2':
                def index_dir = BUILD_BOWTIE2_INDEX(ref_fasta)
                ALIGN_BOWTIE2(ch_trimmed_for_align, index_dir.collect())
                ch_aligned_with_meta = ALIGN_BOWTIE2.out.bam_output 
                ch_flagstats = ALIGN_BOWTIE2.out.flagstat 
                break
            case 'bwa':
                def index_dir = BUILD_BWA_INDEX(ref_fasta)
                ALIGN_BWA(ch_trimmed_for_align, index_dir.collect())
                ch_aligned_with_meta = ALIGN_BWA.out.bam_output 
                ch_flagstats = ALIGN_BWA.out.flagstat 
                break
            default:
                error "Unknown aligner: ${params.aligner}"
        }
    }

    // --- HISTORIC BRIDGE: Unwrap here ---
    // Changed .set to assignment to satisfy DSL2 compiler
    def aligned_bams = ch_aligned_with_meta
        .map { meta, bam, bai -> [ meta.id, bam, bai ] }

    // --- Optional STEP 5: Primer masking ---
    def primer_bed_path = params.primer_file ? file(params.primer_file).toAbsolutePath() : null
    if(params.mask_primers || (params.primer_file && params.mask_primers != false)) {
        MASK_PRIMERS(aligned_bams.combine(Channel.value(primer_bed_path)))
        aligned_bams = MASK_PRIMERS.out.mask_primers_output
    }

    // --- Optional STEP 6 & 7: Identity & SMOR ---
    if(params.identity) {
        def result = IDENTITY_FILTER(aligned_bams)
        aligned_bams = result[0]
    }
    if(params.smor) {
        def result = SMOR(aligned_bams)
        aligned_bams = result[0]
    }

    // --- STEP 8: Forking ---
    // Changed .set to assignment to fix "Missing name" error
    def ch_split = aligned_bams.multiMap { sample_id, bam, bai ->
        asap: [ sample_id, bam, bai ]   
        ivar: [ [id: sample_id], bam, bai ] 
    }
    
    // --- STEP 9: ASAP Processing ---
    if(params.asap_snps) {
        def xml_output = PROCESS_BAM(ch_split.asap.combine(json_ch))
        
        // --- STEP 9.1: ASAP Tools Processing ---
        if(params.asaptools_processing) {
            PROCESS_XML_R(xml_output, params.proportion)
        }

        if(params.combine_output) {
            def xmls = xml_output.map { id, f -> f }.collect()
            def final_xml = OUTPUT_COMBINER(xmls)
            def default_stylesheet = file("${workflow.projectDir}/../output_transforms/ASAP_fulldetails_web.xsl")
            def stylesheet_path = params.stylesheet ? file(params.stylesheet) : default_stylesheet
            FORMAT_OUTPUT(final_xml, Channel.value(stylesheet_path))
        }
    }

    // --- STEP 10: iVAR Trimming ---
    def ch_bam_for_ivar
    if (params.ivar || params.ivar_trim) {
        IVAR_TRIM (ch_split.ivar, primer_bed_path)
        ch_bam_for_ivar = IVAR_TRIM.out.bam
    } else {
        ch_bam_for_ivar = ch_split.ivar.map { meta, bam, bai -> [meta, bam] }
    }

    // --- STEP 11 & 12: iVAR --- 
    if (params.ivar || params.ivar_variants) { IVAR_VARIANTS (ch_bam_for_ivar, ref_fasta, true) }
    if (params.ivar || params.ivar_consensus) { IVAR_CONSENSUS (ch_bam_for_ivar, ref_fasta, true) }

    // --- STEP 13: MultiQC ---
    MULTIQC (
        FASTQC_INITIAL.out.zip.map{ it[1] }.mix(FASTQC_POST.out.zip.map{ it[1] })
            .mix(ch_trim_json_for_multiqc.map{ it[1] })
            .mix(ch_flagstats.map{ it[1] })
            .collect(),
        [], [], [], [], []
    )
}

// Sub-workflows
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