process CONVERT_TO_FASTQ {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir path: "${params.outdir}/fastqs/", enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*.fastq.gz"
    
    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(reads_bam)

    output:
    tuple val(meta), path(fastq_output),  emit: reads_fastq

    script:
    //will likely need thinking here as if we do other methods of sequencing and use this pipeline need the correct out flags
    fastq_output = "${meta.ID}.fastq.gz"
    """
    samtools fastq -@ ${task.cpus} -0 ${fastq_output} ${reads_bam}
    """
}

process MERGE_BAMS_FOR_SUMMARY {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    path("*.bam")

    output:
    path(combined_bam),  emit: summary_bam

    script:
    combined_bam = "merged.bam"
    """
    samtools merge -@ ${task.cpus} -o ${combined_bam} *.bam
    """
}

process MANAGE_DUPLICATES_FROM_BAMS {
    label 'cpu_2'
    label 'mem_4'
    label 'time_1'

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(bam), path(duplicates_list)

    output:
    tuple val(meta), path(final_bam), emit: bam

    script:
    final_bam = "${bam.simpleName}_clean.bam"
    """
    samtools view -N ${duplicates_list} -o ${final_bam} ${bam}
    """
}

process CONVERT_TO_BAM {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_1'
    
    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(mapped_reads)

    output:
    tuple val(meta), path("${mapped_reads_bam}"),  emit: mapped_reads_bam

    script:
    mapped_reads_bam = "${meta.ID}.bam"
    """
    samtools view \
        -@ ${task.cpus} \
        -bS -F 2308 \
        -o ${mapped_reads_bam} \
        ${mapped_reads}
    """
}

process SAMTOOLS_SORT {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/mapped_reads", enabled: params.keep_sorted_bam, mode: 'copy', overwrite: true

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(mapped_reads_bam)

    output:
    tuple val(meta), path("${sorted_reads}"),  emit: sorted_reads

    script:
    sorted_reads = "${meta.ID}_sorted.bam"
    """
    samtools sort \
        -@ ${task.cpus} \
        -o ${sorted_reads} \
        ${mapped_reads_bam}
    """
}

process INDEX_REF {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/sorted_ref", mode: 'copy', overwrite: true

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    path(reference)

    output:
    tuple path(reference), path("${faidx}"),  emit: ref_index

    script:
    faidx = "${reference}.fai"
    """
    samtools faidx "${reference}" > "${faidx}"
    """
}

process SAMTOOLS_INDEX_BAM {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    publishDir "${params.outdir}/mapped_reads", enabled: params.keep_bam_files, mode: 'copy', overwrite: true, pattern: "*.bai"

    input:
    tuple val(meta), path(bam_file)
    
    output:
    tuple val(meta), path(bam_file), path("*.bai"),  emit: bam_index
    
    script:
    """
    samtools index -@ ${task.cpus} *.bam
    """
}


process GET_READLENGTH_DISTRIBUTION {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/qc/readlengths", mode: 'copy', overwrite: true

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(sorted_bam), path(bam_index)

    output:
    tuple val(meta), path("*.read-lengths.tsv"),  emit: readlengths

    script:
    sorted_bam = "${meta.ID}_sorted.bam"
    """
    samtools view -@ ${task.cpus} ${sorted_bam} | \
    awk '{print length(\$10)}' | sort | uniq -c | sort -n -k 2 | awk -v OFS='\t' '{print \$2,\$1}' \
    > ${meta.ID}.read-lengths.tsv 
    """
}

process SAMTOOLS_DEPTH {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    publishDir "${params.outdir}/qc/coverage/samtools_depth", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam_file), path(bam_index)

    output:
    tuple val(meta), path(coverage_report),  emit: samtools_coverage

    script:
    coverage_report = "${meta.ID}_samtools_depth.tsv"
    """
    samtools depth -@ ${task.cpus} -aa *.bam -o ${coverage_report}
    """
}
