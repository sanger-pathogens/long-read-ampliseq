process BASECALL {
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'
    label 'gpu'

    container '/data/pam/installs/images/dorado-0.5.1.simg'
    
    input:
    path(pod5)

    output:
    path("calls.bam"), emit: called_channel

    script:
    basecall_model = "${params.basecall_model == "auto" ? "" : ",${params.basecall_model}"}"
    """
    dorado basecaller ${params.basecalling_accuracy}${basecall_model} --trim ${params.trim_adapters} ${pod5} > calls.bam
    """
}

process DEMUX {
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'
    label 'gpu'

    tag "${barcode_kit_name}"

    container 'quay.io/sangerpathogens/cuda_dorado:0.5.1'
    
    input:
    path(called_bam)
    each(barcode_kit_name)

    output:
    tuple val(barcode_kit_name), path("barcodes/*.bam"), emit: called_channel

    script:
    """
    dorado demux --output-dir ./barcodes --kit-name ${barcode_kit_name} ${called_bam}
    """
}

process DORADO_SUMMARY { 
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'

    publishDir path: "${params.outdir}/sequencing_summary/", enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "summary.tsv"

    container 'quay.io/sangerpathogens/cuda_dorado:0.5.1'
    
    input:
    path(called_bam)

    output:
    path("summary.tsv"), emit: summary_channel

    script:
    """
    dorado summary ${called_bam} > summary.tsv
    """
}