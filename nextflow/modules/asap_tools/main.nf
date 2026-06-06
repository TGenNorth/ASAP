#! /usr/bin/env nextflow

process PROCESS_XML_R {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}/XML_Rdata", mode: 'copy'

    input:
    tuple val(sample_id), path(xml)
    val proportion

    output:
    path "${sample_id}_XML_Data.Rdata", emit: rdata
    path "${sample_id}_Summary.csv",   emit: csv

    script:
    """
    process_xml.R ${xml} ${proportion} ${sample_id}
    """
}

process PROCESS_COMBINE_RDATA {
    tag "combine_rdata"
    label 'process_medium'

    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    val  poi_input
    path rdata_files

    output:
    path "Combined_ASAP_Data.Rdata", emit: combined_rdata
    path "Combined_Summary.csv",   emit: combined_csv

    script:
    def poi_param = (poi_input == null || poi_input == "NULL" || poi_input == "") ? "NULL" : poi_input
    """
    process_combine_rdata.R ${poi_param} ${rdata_files}
    """
}

process PROCESS_GENERATE_FASTA {
    tag "GENERATE_FASTA"
    label 'process_low'

    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    path combined_rdata
    val  prefix

    output:
    path "*.fasta", emit: fasta, optional: true

    script:
    def fasta_threshold = params.asaptools_breadth_threshold ?: params.breadth
    """
    process_asaptools_fasta_export.R ${combined_rdata} ${prefix} ${fasta_threshold}
    """
}

process PROCESS_GENERATE_COV_TABLE {
    tag "coverage_table"
    label 'process_medium'

    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    path combined_rdata
    val  min_depth
    val  prefix
    val  poi_input

    output:
    path "*.xlsx", emit: excel, optional: true

    script:
    def poi_param = (poi_input == null || poi_input == "NULL") ? "NULL" : poi_input
    """
    process_asaptools_cov_table.R ${combined_rdata} ${min_depth} ${prefix} ${poi_param}
    """
}

process PROCESS_SNPS_TO_AMINOACIDS {
    tag "snp_to_aa_conversion"

    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    path combined_rdata
    path "genbank_input/*"

    output:
    path "SNP_Amino_Acid_Table.Rdata", emit: snp_to_amino_rdata
    path "*.csv",  emit: csv,   optional: true

    script:
    """
    process_asaptools_snps_amino_acids.R \\
        ${combined_rdata} \\
        genbank_input/*
    """
}

process PROCESS_GENERATE_SNP_TABLE {
    tag "snp_table"
    label 'process_medium'

    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    path combined_rdata
    val  prefix
    path "genbank_input/*"
    path primer_bed
    val  poi_input
    path aa_rdata

    output:
    path "*.xlsx", emit: xlsx, optional: true
    path "*.csv", emit: csv, optional: true

    script:
    def poi_param = (poi_input == null || poi_input == "NULL" || poi_input == "") ? "NULL" : poi_input
    def exclude_list = (params.asaptools_samples_to_remove == null || params.asaptools_samples_to_remove == "") ? "NONE" : params.asaptools_samples_to_remove
    def effective_prop = params.asaptools_snp_proportion ?: params.proportion
    def xls_toggle = params.asaptools_snp_table_xls.toString().toUpperCase()
    def bed_param = (primer_bed && primer_bed.name != 'null') ? primer_bed : "NULL"
    def aa_param  = (aa_rdata && aa_rdata.name != 'null') ? aa_rdata : "NULL"

    """
    shopt -s nullglob
    process_asaptools_snp_table.R \\
        ${combined_rdata} \\
        ${prefix} \\
        ${effective_prop} \\
        ${params.asaptools_max_sample_snp_count} \\
        ${params.asaptools_min_location_depth} \\
        "${exclude_list}" \\
        ${poi_param} \\
        "${bed_param}" \\
        "${aa_param}" \\
        ${xls_toggle} \\
        genbank_input/*
    """
}

process PROCESS_QC_PLOTS {
    tag "qc_plots"

    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    path combined_rdata
    val  prefix
    val  poi_input

    output:
    path "*.html", emit: html
    path "*.jpg",  emit: jpg

    script:
    def poi_param = (poi_input == null || poi_input == "NULL" || poi_input == "") ? "NULL" : poi_input
    """
    process_asaptools_generate_figures.R \\
        ${combined_rdata} \\
        ${prefix} \\
        ${poi_param} \\
        ${params.asaptools_snp_proportion} \\
        ${params.asaptools_min_location_depth}
    """
}
