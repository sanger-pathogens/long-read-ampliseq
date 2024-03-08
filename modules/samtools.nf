process CONVERT_TO_FASTQ {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), file(reads_bam)

    output:
    tuple val(meta), path("*.fastq.gz"),  emit: reads_fastq

    script:
    //will likely need thinking here as if we do other methods of sequencing and use this pipeline need the correct out flags
    """
    samtools fastq \
        -@ ${task.cpus} \
        -0 ${meta.ID}.fastq.gz \
        -1 ${meta.ID}_1.fastq.gz \
        -2 ${meta.ID}_2.fastq.gz \
        ${reads_bam}
    """
}