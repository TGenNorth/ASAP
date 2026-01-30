process FASTPLONG {
    tag "$meta.id"
    label 'process_medium'
    // Ensure results are saved to the results folder like RUN_FASTP
    publishDir "${params.outdir}/sample_info/${meta.id}/fastp_long", mode: 'copy'

    container "quay.io/biocontainers/fastplong:0.4.1--h224cc79_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastplong.fastq.gz"), emit: reads
    path "*.fastplong.json"                      , emit: json
    path "*.fastplong.html"                      , emit: html
    path "versions.yml"                          , emit: versions

    script:
    // This allows you to use the same parameter you used for short-reads, 
    // or a specialized one if you prefer.
    def extra_args = params.fastp_extra_args ?: ""
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    fastplong \\
        -i $reads \\
        -o ${prefix}.fastplong.fastq.gz \\
        -j ${prefix}.fastplong.json \\
        -h ${prefix}.fastplong.html \\
        --thread $task.cpus \\
        ${extra_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: \$(fastplong --version 2>&1 | sed -e "s/fastplong //g")
    END_VERSIONS
    """
}