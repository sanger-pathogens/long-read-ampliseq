import org.apache.commons.io.FilenameUtils

include { CONVERT_FAST5_TO_POD5 } from '../modules/pod5.nf'
include { BASECALL; DEMUX } from '../modules/dorado.nf'

def validateSingleFormat(listOfFormats){
    if (listOfFormats.size() != 1) {
        log.error("Multiple signal filetypes ${listOfFormats} found in '${params.raw_read_dir}'. Please separate filetypes into distinct directories and process indepedently.")
    }
}

workflow BASECALLING {  
    take:
    raw_read_signal_files

    main:
    raw_read_signal_files.map{ raw_read_signal_file ->
        tuple(raw_read_signal_file.extension, raw_read_signal_file)
    }
    | groupTuple
    | set{ input_formats }

    /* 
    validates that we only have 1 file format in basecall dir otherwise exit
    */

    input_formats.map{ format, files -> format }
    | collect
    | map{ format -> validateSingleFormat(format)}

    /*
    Files in the fast5 format are converted to pod5 and so are branched out into their respective channels
    */

    input_formats
    | branch{format, files ->
        fast5: format == "fast5"
            return files

        pod5: format == "pod5"
           return files
    }
    | set{ raw_files }

    CONVERT_FAST5_TO_POD5(raw_files.fast5)
    | mix(raw_files.pod5) //files that were already pod5 are added back in after the convert process
    | BASECALL
    
    DEMUX(BASECALL.out.called_channel, params.barcode_kit_name) //todo https://github.com/nanoporetech/dorado/issues/625 if list of barcodes provided loop over and call for each
    | transpose
    | map{ barcode_kit, bam -> 
        def meta = [:]
        meta.barcode_kit = barcode_kit
        meta.barcode = "${ bam.simpleName.contains("barcode") ? bam.simpleName.split("barcode")[1] : bam.simpleName }" //i.e. when bam.simpleName = unclassified
        tuple(meta, bam)
    }
    | set{ barcode_bam_ch }

    barcode_bam_ch.view()
    emit:
    barcode_bam_ch

}