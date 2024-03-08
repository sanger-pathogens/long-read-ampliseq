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
    if (params.basecall_model == "")
        """
        dorado basecaller ${params.basecalling_accuracy} ${pod5} --trim ${params.trim_adapters} > calls.bam
        """
    else
        """
        dorado basecaller ${params.basecalling_accuracy},${params.basecall_model} ${pod5} --trim ${params.trim_adapters} > calls.bam
        """
}

process DEMUX {
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'
    label 'gpu'

    container '/data/pam/installs/images/dorado-0.5.1.simg'
    
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