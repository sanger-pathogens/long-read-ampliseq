process NANOPLOT_QC {
    container "quay.io/biocontainers/nanoplot:1.38.0--pyhdfd78af_0"
    publishDir "${params.results_dir}/qc/nanoplot", mode: 'copy', overwrite: true
    input:
        path(fastq_file)
    output:
        path("*_nanoplot_qc/*")
    script:
        """
        sample=\$(basename $fastq_file | awk -F "." '{ print \$1}')
        NanoPlot -t 2 --fastq $fastq_file -o \${sample}_nanoplot_qc
        """
}