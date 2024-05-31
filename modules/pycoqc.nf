process PYCOQC {
    container "quay.io/biocontainers/pycoqc:2.5.2--py_0"
    publishDir "${params.outdir}/qc/pycoqc", mode: 'copy', overwrite: true
    
    input:
        path(sequence_summary_file)

    output:
        path("*.html"), emit: html
        path("*.json"), emit: json
        
    script:
    final_qc_file = "summary_pycoqc"
    """
    pycoQC -f ${sequence_summary_file} -o ${final_qc_file}.html -j ${final_qc_file}.json
    """
}