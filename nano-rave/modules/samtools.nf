process GET_CHROM_SIZES_AND_INDEX {
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    input:
        tuple val(ref_id), path(reference)
    output:
        tuple val(ref_id), path("*.fai"), emit: fasta_index_ch
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        samtools faidx ${reference}
        """
}

process SAMTOOLS_VIEW_SAM_TO_BAM {
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    input:
        tuple val(ref_id), path(reference)
        path(sam_file)
    output:
        path("*.bam"), emit: bam_file
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        filename=\$(basename $sam_file | awk -F "." '{ print \$1}')
        samtools view -b -h -O BAM -@ 2 -o \${filename}.bam $sam_file
        """
}

process SAMTOOLS_SORT_AND_INDEX {
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    if (params.keep_bam_files) {
        publishDir "${params.results_dir}/bams", mode: 'copy', overwrite: true, pattern: "*.bam*"
    }
    input:
        tuple val(ref_id), path(reference)
        path(bam_file)
    output:
        tuple path("*.sorted.bam"), path("*.sorted.bam.bai"),  emit: sorted_bam
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        filename=\$(basename $bam_file | awk -F "." '{ print \$1}')
        samtools sort -@ 2 -o \${filename}.sorted.bam -T \${filename} $bam_file
        samtools index \${filename}.sorted.bam
        """
}