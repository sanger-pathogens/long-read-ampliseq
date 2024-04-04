process GUNZIP {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'

    input:
    tuple val(meta), path(zipped_gvcf), path(reference), path(reference_index)

    output:
    tuple val(meta), path("merge_output.gvcf"), path(reference), path(reference_index), emit: gunzip

    script:
    """
    gunzip -c ${zipped_gvcf} > merge_output.gvcf
    """
}