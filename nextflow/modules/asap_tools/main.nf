#! /usr/bin/env nextflow

process PROCESS_XML_R {
    tag "$sample_id"
    label 'process_low'
    
    // Direct path to your existing environment
    conda "/tgen_labs/EPIC/miniconda3/envs/r_mirror_env"

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
    conda "/tgen_labs/EPIC/miniconda3/envs/r_mirror_env"
    
    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    val  poi_input
    path rdata_files // The list of all .Rdata files from .collect()

    output:
    path "Combined_ASAP_Data.Rdata", emit: combined_rdata
    path "Combined_Summary.csv",   emit: combined_csv

    script:
    def poi_param = (poi_input == null || poi_input == "NULL" || poi_input == "") ? "NULL" : poi_input
    """
    # Use a shell script wrapper or call R directly
    process_combine_rdata.R ${poi_param} ${rdata_files}
    """
}

process PROCCESS_GENERATE_FASTA {
    tag "GENERATE_FASTA"
    label 'process_low'
    conda "/tgen_labs/EPIC/miniconda3/envs/r_mirror_env"
    
    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    path combined_rdata  // arg[1]
    val  prefix          // arg[2]

    output:
    path "*.fasta", emit: fasta, optional: true

    script:
    def fasta_threshold = params.asaptools_breadth_threshold ?: params.breadth
    """
    process_asaptools_fasta_export.R ${combined_rdata} ${prefix} ${fasta_threshold}
    """
}


process PROCCESS_GENERATE_COV_TABLE {
    tag "coverage_table"
    label 'process_medium'
    conda "/tgen_labs/EPIC/miniconda3/envs/r_mirror_env"
    
    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    path combined_rdata  // arg[1]
    val  min_depth       // arg[2]
    val  prefix          // arg[3]
    val  poi_input       // arg[4] - Changed from 'any' to 'val'

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
    label 'process_medium'
    conda "/tgen_labs/EPIC/miniconda3/envs/r_mirror_env"

    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    path combined_rdata      // 1
    path genbank_ref        // 3

    output:
    path "SNP_Amino_Acid_Table.Rdata", emit: snp_to_amino_rdata
    path "*.csv",  emit: csv,   optional: true

    script:
    """
    process_asaptools_snps_amino_acids.R \\
        ${combined_rdata} \\
        ${genbank_ref}
    """
}

process PROCESS_GENERATE_SNP_TABLE {
    tag "snp_table"
    label 'process_medium'
    conda "/tgen_labs/EPIC/miniconda3/envs/r_mirror_env"

    publishDir "${params.outdir}/ASAP_R_Data", mode: 'copy'

    input:
    path combined_rdata      // 1
    val  prefix              // 2
    path genbank_ref        // 3
    path primer_bed         // 4
    val  poi_input           // 5
    path aa_rdata           // 6 - This maps to the file from snp_amino_data
    
    output:
    path "*.xlsx", emit: excel, optional: true

    script:
    def poi_param = (poi_input == null || poi_input == "NULL" || poi_input == "") ? "NULL" : poi_input
    def exclude_list = (params.asaptools_samples_to_remove == null || params.asaptools_samples_to_remove == "") ? "NONE" : params.asaptools_samples_to_remove
    
    """
    process_asaptools_snp_table.R \\
        ${combined_rdata} \\
        ${prefix} \\
        ${params.asaptools_snp_proportion} \\
        ${params.asaptools_max_sample_snp_count} \\
        ${params.asaptools_min_location_depth} \\
        "${exclude_list}" \\
        ${poi_param} \\
        ${genbank_ref} \\
        ${primer_bed} \\
        ${aa_rdata} 
    """
}