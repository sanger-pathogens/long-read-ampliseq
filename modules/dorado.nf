process MODEL_DOWNLOAD {
    label 'cpu_1'
    label 'mem_4'
    label 'time_30m'

    container 'quay.io/sangerpathogens/cuda_dorado:0.5.1'
    
    input:
    path(pod5)

    output:
    tuple path(pod5), path(basecall_model), emit: model_ch

    script:
    basecall_model = "${params.basecall_model}"
    """
    dorado download --model ${basecall_model}
    """
}

process BASECALL {
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'
    label 'gpu'

    container 'quay.io/sangerpathogens/cuda_dorado:0.5.1'
    
    input:
    tuple path(pod5), path(model)

    output:
    path("calls.bam"), emit: called_channel

    script:
    min_qscore = "${params.min_qscore == "" ? "" : "--min-qscore ${params.min_qscore}"}"
    """
    dorado basecaller ${model} --trim ${params.trim_adapters} ${min_qscore} ${pod5} > calls.bam
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

    publishDir path: "${params.outdir}/sequencing_summary/", mode: 'copy', overwrite: true, pattern: "summary.tsv"

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

process UNASSIGNED_SUMMARY { 
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'

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