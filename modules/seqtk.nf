process MASK_READS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/biocontainers/seqtk:1.4--he4a0461_2'

    input:
    tuple val(meta), path(fastq_trimmed)

    output:
    tuple val(meta), path("${fastq_trimmed_split}.gz"), emit: trimmed_split_reads

    script:
    fastq_trimmed_split = "${meta.ID}_trimmed_split.fastq"
    """
    seqtk seq -q 15 -n N ${fastq_trimmed} > ${fastq_trimmed_split}

    gzip ${fastq_trimmed_split}
    """
}
