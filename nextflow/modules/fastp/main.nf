#! /usr/bin/env nextflow

process RUN_FASTP {
    tag "${meta.id}"
    stageInMode 'copy'
    publishDir "${params.outdir}/sample_info/${meta.id}/fastp", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path adapter_fasta

    output:
    // FIX: The asterisk before AND after 'cleaned' ensures SE and PE are both caught
    tuple val(meta), path("*.cleaned*.fastq.gz"), emit: trimmed_reads
    path "*.fastp.html",                        emit: html
    path "*.fastp.json",                        emit: json

    script:
    def extra_args = params.fastp_extra_args ?: ""
    // Define names explicitly to avoid globbing confusion
    def out1 = "${meta.id}.cleaned_R1.fastq.gz"
    def out2 = "${meta.id}.cleaned_R2.fastq.gz"
    def out_se = "${meta.id}.cleaned.fastq.gz"

    if (meta.single_end) {
        """
        fastp -i ${reads[0]} -o ${out_se} --adapter_fasta ${adapter_fasta} \\
            --html ${meta.id}.fastp.html --json ${meta.id}.fastp.json --thread ${task.cpus} ${extra_args}
        """
    } else {
        """
        fastp -i ${reads[0]} -I ${reads[1]} -o ${out1} -O ${out2} --adapter_fasta ${adapter_fasta} \\
            --html ${meta.id}.fastp.html --json ${meta.id}.fastp.json --thread ${task.cpus} ${extra_args}
        """
    }
}