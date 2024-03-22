process CUT_PRIMERS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir path: "${params.outdir}/cutadapt/", enabled: params.save_trimmed, mode: 'copy', overwrite: true, pattern: "*_trimmed.fastq.gz"
    publishDir path: "${params.outdir}/cutadapt/", enabled: params.save_too_short, mode: 'copy', overwrite: true, pattern: "*_too_short.fastq.gz"
    publishDir path: "${params.outdir}/cutadapt/", enabled: params.save_too_long, mode: 'copy', overwrite: true, pattern: "*_too_long.fastq.gz"
    
    conda 'bioconda::cutadapt=4.7'
    container 'quay.io/biocontainers/cutadapt:4.7--py310h4b81fae_1'

    input:
    tuple val(meta), path(fastq), path(primers)

    output:
    tuple val(meta), path(fastq_trimmed), emit: trimmed_reads
    tuple val(meta), path(too_short_reads), emit: too_short_reads
    tuple val(meta), path(too_long_reads), emit: too_long_reads

    script:
    fastq_trimmed = "${meta.ID}_trimmed.fastq.gz"
    too_short_reads = "${meta.ID}_too_short.fastq.gz"
    too_long_reads = "${meta.ID}_too_long.fastq.gz"
    // TODO: As noted in Mat's script, -g looks at 5', whilst -a looks at 3' and -b looks at both (but seems overly aggressive).
    // Could look at trimming using 'linked primers', but would need separate primer files for each end.
    """
    cutadapt \
        -j ${task.cpus} \
        -g "file:${primers}" \
        -o ${fastq_trimmed} \
        -m ${params.lower_read_length_cutoff} \
        -M ${params.upper_read_length_cutoff} \
        --too-short-output ${too_short_reads} \
        --too-long-output ${too_long_reads} \
        ${params.cutadapt_args} \
        *.fastq.gz
    """
}

