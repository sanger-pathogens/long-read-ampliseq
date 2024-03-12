process PYCOQC {
    container "quay.io/biocontainers/pycoqc:2.5.2--py_0"
    publishDir "${params.results_dir}/qc/pycoqc", mode: 'copy', overwrite: true
    input:
        tuple val(sequencing_dir), val(sequence_summary_file)
    output:
        path("*.html")
        path("*.json")
    script:
        """
        sample=\$(echo $sequencing_dir | awk -F "/" '{ print \$(NF-1) }')
        pycoQC -f $sequence_summary_file -o \${sample}_pycoqc.html -j \${sample}_pycoqc.json
        """
}