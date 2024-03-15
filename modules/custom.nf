process PYTHON_COVERAGE_OVER_DEFINED_REGIONS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    conda "bioconda::pandas=2.2.1"
    container "quay.io/sangerpathogens/pandas:2.2.1"
    // conda "bioconda::pandas=1.5.2"
    // container "quay.io/biocontainers/pandas:1.5.2"

    publishDir "${params.outdir}/${meta.ID}/qc/samtools_coverage", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam_file), path(bam_index), path(samtools_coverage), path(target_regions_bed)
    
    output:
    tuple val(meta), path(coverage_summary),  emit: coverage_summary
    
    script:
    coverage_summary = "${meta.ID}_sorted.depth.tsv"
    """
    coverage_over_defined_regions.py \
        -s ${samtools_coverage} \
        -b ${bam_file} \
        -c ${target_regions_bed} \
        -w
    """
}