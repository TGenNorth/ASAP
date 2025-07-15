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
    bwa index -p bwa_index ${reference_fasta}
    mkdir bwa_index
    cp bwa_index.* bwa_index/
    """
}

// Align read pairs using BWA index
process ALIGN_BWA {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}/bwa", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2), path(index_dir)

    output:
    tuple val(sample_id), path("${sample_id}-bwa.bam"), path("${sample_id}-bwa.bam.bai"), emit: bam_output

    script:
    """
    # Copy all the index files to local working dir
    cp ${index_dir}/* .

    # Run bwa mem on the paired reads
    bwa mem -t ${task.cpus} -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' bwa_index ${read1} ${read2} \
    | samtools view -Sbh - \
    | samtools sort -T ${sample_id}-bwa -o ${sample_id}-bwa.bam - 

    # Index the bam
    samtools index ${sample_id}-bwa.bam
    """
}

