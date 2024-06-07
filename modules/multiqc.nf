process MULTIQC {
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    //TODO Check if user changes to supplied multiqc config are recognised without this. They didn't seem to be when I was editing the pipeline default multiqc_config.yml.
    cache false

    conda 'bioconda::multiqc=1.21'
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    path('pycoqc/*')
    path('fastqc/*')
    path('samtools_stats/*')

    output:
    path("multiqc_report.html"), emit: report
    path("*_data"), emit: data
    path("*_plots"), optional:true, emit: plots

    script:
    def custom_config = params.multiqc_config ? "--config ${params.multiqc_config}" : "--config ${projectDir}/config/multiqc/multiqc_config.yml"
    """
    multiqc \
        -n multiqc_report.html \
        -f \
        ${custom_config} \
        .
    """
}