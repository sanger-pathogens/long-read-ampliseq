process MERGE_GVCF {
    label 'cpu_1'
    label 'mem_2'
    label 'time_12'

    publishDir "${params.outdir}/variants/merged_gvcf/", mode: 'copy', overwrite: true, pattern: "*_merged.vcf.gz"

    conda "bioconda::bcftools=1.20"
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"

    input:
    path("variant_vcfs*.gvcf.gz")

    output:
    path("${workflow.runName}_merged.vcf.gz"), emit: merged_vcf

    script:
    //workflow start is ugly 2024-02-29T12:01:26.233465Z split on T for time and take just the date
    date = "${workflow.start}".split('T')[0]
    """
    find . -name "*.gvcf.gz" -exec bcftools index "{}" \\;

    bcftools merge -R ${params.target_regions_bed} -0 -o ${date}_merged.vcf.gz -O z -f PASS *.gvcf.gz
    """
}