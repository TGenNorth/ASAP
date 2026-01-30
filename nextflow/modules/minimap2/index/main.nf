process MINIMAP2_INDEX {
    tag "$meta.id"
    label 'process_low'

    container "community.wave.seqera.io/library/minimap2_samtools:f2e55a1cd407fcaf"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -d ${fasta.baseName}.mmi \\
        -t $task.cpus \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}