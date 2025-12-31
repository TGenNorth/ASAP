#! /usr/bin/env nextflow

process RUN_FASTP {
    tag "$sample_id"
    stageInMode 'copy'
    publishDir "${params.outdir}/${sample_id}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path adapter_fasta

    output:
    // This emits exactly 5 items: [val, path, path, path, path]
    tuple val(sample_id), path("${sample_id}.cleaned_R1.fastq.gz"), path("${sample_id}.cleaned_R2.fastq.gz"), path("${sample_id}.fastp.html"), path("${sample_id}.fastp.json"), emit: trimmed_reads
    path "${sample_id}.fastp.html", emit: html
    path "${sample_id}.fastp.json", emit: json

    script:
    def extra_args = params.fastp_extra_args ?: ""
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${sample_id}.cleaned_R1.fastq.gz \\
        -O ${sample_id}.cleaned_R2.fastq.gz \\
        --adapter_fasta ${adapter_fasta} \\
        --html ${sample_id}.fastp.html \\
        --json ${sample_id}.fastp.json \\
        --thread ${task.cpus} \\
        ${extra_args}
    """
}

