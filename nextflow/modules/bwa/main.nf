#! /usr/bin/env nextflow

// Build BWA index from reference.fasta
process BUILD_BWA_INDEX {
    tag "bwa_index"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path reference_fasta

    output:
    path 'bwa_index', emit: index_dir

    script:
    """
    mkdir bwa_index
    bwa index -p bwa_index/genome ${reference_fasta}
    """
}

// Align reads using BWA index
process ALIGN_BWA {
    tag "${meta.id}"
    publishDir "${params.outdir}/sample_info/${meta.id}/bwa", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index_dir

    output:
    tuple val(meta), path("${meta.id}-bwa.bam"), path("${meta.id}-bwa.bam.bai"), emit: bam_output
    path "${meta.id}.flagstat.txt", emit: flagstat

    script:
    def extra_args = params.aligner_extra_args ?: ""
    // In Groovy, a list (reads) joined by a space handles one or two files perfectly for BWA
    def input_reads = reads.join(' ')
    
    """
    # Run bwa mem
    # Note: BWA handles both SE and PE by just passing the file(s) at the end
    bwa mem -t ${task.cpus} \\
        -R '@RG\\tID:${meta.id}\\tSM:${meta.id}' \\
        ${extra_args} \\
        ${index_dir}/genome \\
        ${input_reads} \\
    | samtools view -Sbh - \\
    | samtools sort -T ${meta.id}-bwa -o ${meta.id}-bwa.bam - 

    # Index the bam
    samtools index ${meta.id}-bwa.bam
    
    # Generate stats for MultiQC
    samtools flagstat ${meta.id}-bwa.bam > ${meta.id}.flagstat.txt
    """
}