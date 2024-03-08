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
        tuple(FilenameUtils.getExtension(raw_read_signal_file.toString()), raw_read_signal_file)
    }
    | groupTuple
    | set{ input_formats }

    /* 
    validates that we only have 1 file format in basecall dir otherwise exit
    */

    input_formats.map{ format, files -> format }
    | collect
    | map{ format -> validateSingleFormat(format)}

    input_formats.branch{format, files ->
        fast5: format == "fast5"
            return files

        pod5: format == "pod5"
           return files
    }.set{ raw_files }

    CONVERT_FAST5_TO_POD5(raw_files.fast5)
    | mix(raw_files.pod5)
    | BASECALL
    | DEMUX
    | flatten
    | map{ bam -> tuple(FilenameUtils.getBaseName(bam.name), bam)}
    | set{ barcode_bam_ch }

    emit:
    barcode_bam_ch

}