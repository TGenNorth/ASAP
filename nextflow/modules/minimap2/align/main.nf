process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // Nextflow will automatically convert this to a Singularity image
    container "community.wave.seqera.io/library/minimap2_samtools:f2e55a1cd407fcaf"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    val bam_format
    val bam_index_extension
    val cigar_paf_format
    val cigar_bam

    output:
    // Matches BWA exactly: [meta, bam, bai]
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam_output
    path "versions.yml"                               , emit: versions

    script:
    def args   = params.aligner_extra_args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    minimap2 \\
        $args \\
        --MD \\
        -t $task.cpus \\
        $reference \\
        $reads \\
        -a | samtools sort -@ ${task.cpus-1} -o ${prefix}.bam -
    
    samtools index ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
