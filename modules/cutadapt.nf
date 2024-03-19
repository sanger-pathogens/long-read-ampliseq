process CUT_PRIMERS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir path: "${params.outdir}/cutadapt/", enabled: params.save_trimmed, mode: 'copy', overwrite: true, pattern: "*_trimmed.fastq.gz"
    
    conda 'bioconda::cutadapt=4.7'
    container 'quay.io/biocontainers/cutadapt:4.7--py310h4b81fae_1'

    input:
    tuple val(meta), path(fastq), path(primers)

    output:
    tuple val(meta), path(fastq_trimmed),  emit: trimmed_reads

    script:
    fastq_trimmed = "${meta.ID}_trimmed.fastq.gz"
    // TODO: Filter empty reads with -m <minimum_length>? Or save that for a distinct step later...
    // TODO: As noted in Mat's script, -g looks at 5', whilst -a looks at 3' and -b looks at both (but seems overly aggressive).
    // Could look at trimming using 'linked primers', but would need separate primer files for each end.
    // Need to decide on options!
    """
    cutadapt \
        -j ${task.cpus} \
        -g "file:${primers}" \
        -o ${fastq_trimmed} \
        ${params.cutadapt_args} \
        *.fastq.gz
    """
}

