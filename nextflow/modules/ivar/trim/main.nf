process IVAR_TRIM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ivar:1.4.4--h077b44d_0' :
        'biocontainers/ivar:1.4.4--h077b44d_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("${prefix}.bam")    , emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: bai
    tuple val(meta), path("${prefix}.ivar.log"), emit: log
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Define prefix here so it can be used in the output block
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ivar trim \\
        -i $bam \\
        -b $bed \\
        -p ${prefix}_temp \\
        $args \\
        > ${prefix}.ivar.log
        
    # Sort the ivar output (which is currently unsorted)
    samtools sort \\
        -o ${prefix}.bam \\
        ${prefix}_temp.bam

    # samtools index creates ${prefix}.bam.bai
    samtools index ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version | sed -n 's|iVar version \\(.*\\)|\\1|p')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ivar.log
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version | sed -n 's|iVar version \\(.*\\)|\\1|p')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}