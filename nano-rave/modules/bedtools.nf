process BEDTOOLS_GENOMECOV {
    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    publishDir "${params.results_dir}/genome_coverage", mode: 'copy', overwrite: true, pattern: "*.bedGraph"
    input:
        tuple path(sorted_bam_file), path(sorted_bam_index)

    output:
        path("*.bedGraph"), emit: genome_cov_ch

    script:
        """
        filename=\$(basename $sorted_bam_file | awk -F "." '{ print \$1}')
        bedtools genomecov -split -ibam $sorted_bam_file -bga | bedtools sort > \${filename}.bedGraph
        """
}