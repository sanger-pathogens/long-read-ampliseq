process MINIMAP2_INDEX {
    label 'cpu_2'
    label 'mem_4'
    label 'time_1'

    conda "bioconda::minimap2=2.26"
    container "quay.io/biocontainers/minimap2:2.26--he4a0461_2"

    input:
    path(reference)

    output:
    tuple path(reference), path("*.mmi"), emit: ref_mm2_index
        
    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} -d ${reference}.mmi ${reference}
    """
}

process MINIMAP2_ALIGN {
    label 'cpu_2'
    label 'mem_4'
    label 'time_1'

    conda "bioconda::minimap2=2.26"
    container "quay.io/biocontainers/minimap2:2.26--he4a0461_2"

    input:
    tuple val(meta), path(fastq), path(reference), path(mm2_index)

    output:
    tuple val(meta), path("*.sam"), emit: mapped_reads_sam

    script:
    """
    minimap2 \
        -ax map-ont \
        -t ${task.cpus} \
        -a ${mm2_index} \
        ${fastq} \
        > ${meta.ID}.sam
    """
}
