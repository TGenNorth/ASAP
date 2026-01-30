#! /usr/bin/env nextflow

process BUILD_BOWTIE2_INDEX {
    tag "bowtie2_index"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path reference_fasta

    output:
    path 'bt2_index', emit: index_dir

    script:
    """
    mkdir bt2_index
    bowtie2-build ${reference_fasta} bt2_index/bt2_index
    """
}

process ALIGN_BOWTIE2 {
    tag "${meta.id}"
    publishDir "${params.outdir}/sample_info/${meta.id}/bowtie2", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index_dir

    output:
    tuple val(meta), path("${meta.id}-bt2.bam"), path("${meta.id}-bt2.bam.bai"), emit: bam_output
    path "${meta.id}.flagstat.txt", emit: flagstat

    script:
    def extra_args = params.aligner_extra_args ?: ""
    // Determine if we use -U (unpaired) or -1/-2 (paired)
    def input_reads = meta.single_end ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    
    """
    # Run bowtie2 using conditional input string
    bowtie2 -x ${index_dir}/bt2_index \\
        --rg-id '${meta.id}' --rg 'SM:${meta.id}' \\
        ${input_reads} \\
        -p ${task.cpus} \\
        ${extra_args} \\
    | samtools view -Sbh - \\
    | samtools sort -T ${meta.id}-bt2 -o ${meta.id}-bt2.bam -

    # Index the bam file
    samtools index ${meta.id}-bt2.bam

    # Generate stats for MultiQC
    samtools flagstat ${meta.id}-bt2.bam > ${meta.id}.flagstat.txt
    """
}