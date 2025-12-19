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
    bowtie2-build ${reference_fasta} bt2_index
    mkdir bt2_index
    mv bt2_index*.bt2 bt2_index/
    """
}

process ALIGN_BOWTIE2 {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bowtie2", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2), path(index_dir)

    output:
    tuple val(sample_id), path("${sample_id}-bt2.bam"), path("${sample_id}-bt2.bam.bai"), emit: bam_output
    path "${sample_id}.flagstat.txt", emit: flagstat

    script:
    def extra_args = params.aligner_extra_args ? params.aligner_extra_args : ""

    """
    # Copy all index files to local working dir
    cp ${index_dir}/* .
    
    # Run bowtie2 on the paired reads
    bowtie2 -x bt2_index --rg-id '${sample_id}' --rg 'SM:${sample_id}' -1 ${read1} -2 ${read2} -p ${task.cpus} ${extra_args} \\
    | samtools view -Sbh - \\
    | samtools sort -T ${sample_id}-bt2 -o ${sample_id}-bt2.bam -

    # Index the bam file
    samtools index ${sample_id}-bt2.bam

    # Generate stats for MultiQC
    samtools flagstat ${sample_id}-bt2.bam > ${sample_id}.flagstat.txt
    """
}