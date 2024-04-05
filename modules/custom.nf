process PYTHON_COVERAGE_OVER_DEFINED_REGIONS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    conda "bioconda::pandas=2.2.1"
    container "quay.io/sangerpathogens/pandas:2.2.1"
    // conda "bioconda::pandas=1.5.2"
    // container "quay.io/biocontainers/pandas:1.5.2"

    publishDir "${params.outdir}/qc/coverage/coverage_summary", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(samtools_coverage), path(target_regions_bed)
    
    output:
    tuple val(meta), path(coverage_summary),  emit: coverage_summary
    
    script:
    coverage_summary = "*coverage_summary.tsv"
    """
    coverage_over_defined_regions.py \
        -s ${samtools_coverage} \
        -b ${target_regions_bed} \
        -t ${params.coverage_reporting_thresholds},${params.coverage_filtering_threshold} \
        -n ${meta.ID}
    """
}