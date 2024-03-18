process FASTQC {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_2'
    label 'time_12'

    publishDir "${params.outdir}/qc/fastqc", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip

    script:
    """
    fastqc \
        -f fastq \
        --threads ${task.cpus} \
        ${reads}
    """
}
