#! /usr/bin/env nextflow

process PROCESS_XML_R {
    tag "$sample_id"
    label 'process_low'
    
    // Direct path to your existing environment
    conda "/tgen_labs/EPIC/miniconda3/envs/asap_r_env"

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