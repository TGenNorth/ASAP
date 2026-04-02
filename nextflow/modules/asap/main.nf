#! /usr/bin/env nextflow

process GENERATE_REFERENCE_FASTA {
    tag "generate_reference"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path assay_json

    output:
    path "reference.fasta", emit: ref_fasta

    script:
    """
    assayInfo.py ${assay_json} > reference.fasta
    """
}

process MASK_PRIMERS {
    tag "$sample_id"
    publishDir "${params.outdir}/sample_info/${sample_id}/mask_primers", mode: 'copy'

    input:
    tuple val(sample_id), path(bamfile), path(bamindex), path(primer_file)

    output:
    tuple val(sample_id), path("${bamfile.getBaseName()}_primerMasked.bam"), path("${bamfile.getBaseName()}_primerMasked.bam.bai"), emit: mask_primers_output
    tuple val(sample_id), path("primer_masking.tsv"), path("primer_masking.log"), emit: mask_primers_logging

    script:
    def mask_bam_string = params.mask_bam ? "--mask-bam" : "--no-mask-bam"
    def ponly_string = params.primer_only ? "--primer-only" : "--no-primer-only"
    """
    maskPrimers.py -b ${bamfile} -p ${primer_file} --wiggle ${params.wiggle} ${mask_bam_string} ${ponly_string} 
    """
}

process IDENTITY_FILTER {
    tag "$sample_id"
    publishDir "${params.outdir}/sample_info/${sample_id}/identity_filter", mode: 'copy'

    input:
    tuple val(sample_id), path(bamfile), path(bamindex)
    
    output:
    tuple val(sample_id), path("${bamfile.getBaseName()}_identityFiltered.bam"), path("${bamfile.getBaseName()}_identityFiltered.bam.bai"), emit: identity_filter_output
    tuple val(sample_id), path("identity_filtering.log"), emit: identity_filter_logging

    script:
    """
    identityFilter.py -b ${bamfile} -i ${params.identity}
    """
}

process SMOR {
    tag "$sample_id"
    publishDir "${params.outdir}/sample_info/${sample_id}/smor", mode: 'copy'

    input:
    tuple val(sample_id), path(bamfile), path(bamindex)
    
    output:
    tuple val(sample_id), path("${bamfile.getBaseName()}_SMOR.bam"), path("${bamfile.getBaseName()}_SMOR.bam.bai"), emit: smor_output
    tuple val(sample_id), path("smor_processing.log"), emit: smor_logging

    script:
    """
    generateSMORbam.py -b ${bamfile} -c ${params.fill_character} 
    """
}

process SMOR_CORRECTION {
    tag "$sample_id"
    publishDir "${params.outdir}/sample_info/${sample_id}/smor_correction/", mode: 'copy'

    input:
    tuple val(sample_id), path(bamfile), path(bamindex)
    
    output:
    tuple val(sample_id), path("${bamfile.getBaseName()}_SMOR.bam"), path("${bamfile.getBaseName()}_SMOR.bam.bai"), emit: smor_output
    tuple val(sample_id), path("smor_processing.log"), emit: smor_logging

    script:
    """
    generateSMORbam_correction.py -b ${bamfile} -c ${params.fill_character} 
    """
}

process PROCESS_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/xml", mode: 'copy'

    input:
    tuple val(sample_id), path(bamfile), path(bamindex), path(assay_json)

    output:
    tuple val(sample_id), path("${sample_id}.xml"), emit: xml_output

    script:
    // Create a string variable that is either the flag or empty
    def wg_flag = params.whole_genome ? "--whole-genome" : ""
    
    """
    newBamProcessor.py \\
        -j ${assay_json} \\
        -b ${bamfile} \\
        -d ${params.depth} \\
        --breadth ${params.breadth} \\
        -p ${params.proportion} \\
        -m ${params.mutation_depth} \\
        --min-base-qual ${params.min_base_qual} \\
        --consensus-proportion ${params.consensus_proportion} \\
        --fill-gaps ${params.fill_gaps} \\
        --mark-deletions ${params.mark_deletions} \\
        ${wg_flag} \\
        -o ${sample_id}.xml
    """
}

process OUTPUT_COMBINER {
    tag "output_combiner"
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path xml_files
    
    output:
    path("${params.file_name}_analysis.xml"), emit: final_xml

    script:
    """
    outputCombiner.py -x . -n ${params.file_name}
    """
}

process FORMAT_OUTPUT {
    tag "format_output"
    publishDir "${params.outdir}/", mode: 'copy'

    def out_file = params.out_file ? params.out_file : "ASAP_Report_${params.file_name}.html"

    input:
    path final_xml
    path stylesheet
    
    output:
    path("*.html"), emit: asap_output
    path("${params.file_name}/"), emit : extra_output

    script:
    """
    formatOutput.py -x ${final_xml} -s ${stylesheet} -o ${out_file} 
    """
}
