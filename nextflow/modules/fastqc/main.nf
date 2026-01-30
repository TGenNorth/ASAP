process FASTQC {
    tag "${meta.id}-${meta.status ?: ''}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    publishDir "${params.outdir}/sample_info/${meta.id}/fastqc/${meta.status ?: ''}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip

    script:
    def args   = task.ext.args ?: ''
    // Use meta.status (initial/post_process) in the prefix to avoid filename collisions
    def prefix = "${meta.id}${meta.status ? '.' + meta.status : ''}"
    
    // Logic to handle renaming to avoid collisions and simplify MultiQC tracking
    def old_new_pairs = (reads instanceof Path || reads.size() == 1) ? 
        [[ reads, "${prefix}.fastq.gz" ]] : 
        reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.fastq.gz" ] }
    
    def rename_to     = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ _old_name, new_name -> new_name }.join(' ')

    def memory_in_mb = task.memory ? task.memory.toUnit('MB') / task.cpus : 2000
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)

    """
    printf "%s %s\\n" ${rename_to} | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        --memory ${fastqc_memory} \\
        ${renamed_files}
    """
}