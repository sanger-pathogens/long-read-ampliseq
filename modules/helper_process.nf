process GUNZIP {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'

    input:
    tuple val(meta), path(zipped_vcf), path(reference), path(reference_index)

    output:
    tuple val(meta), path("*.vcf"), path(reference), path(reference_index), emit: gunzip

    script:
    """
    gunzip -c ${zipped_vcf} > ${meta.ID}.vcf
    """
}