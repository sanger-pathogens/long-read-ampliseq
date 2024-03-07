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
        dorado basecaller ${params.basecalling_accuracy} ${pod5} --trim ${params.trim_adapters} --kit-name ${params.barcode_kit_name} > calls.bam
        """
    else
        """
        dorado basecaller ${params.basecalling_accuracy},${params.basecall_model} ${pod5} --trim ${params.trim_adapters} --kit-name ${params.barcode_kit_name} > calls.bam
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

    output:
    path("barcodes/*.bam"), emit: called_channel

    script:
    """
    dorado demux --output-dir ./barcodes --kit-name ${params.barcode_kit_name} ${called_bam}
    """
}