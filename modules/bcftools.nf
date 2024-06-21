process MERGE_GVCF {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_12'

    publishDir "${params.outdir}/variants/merged_gvcf/${workflow.runName}_merged.vcf", mode: 'copy', overwrite: true, pattern: "*_merged.vcf"

    conda "bioconda::bcftools=1.20"
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"

    input:
    path(vcfs)

    output:
    path("${workflow.runName}_merged.vcf"), emit: merged_vcf

    script:
    """
    bcftools merge -R ${target_regions_bed} -0 -o ${workflow.runName}_merged.vcf -O v -f PASS *_clair3.gvcf.gz
    """
}