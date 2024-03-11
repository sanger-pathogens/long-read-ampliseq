import org.apache.commons.io.FilenameUtils

include { CONVERT_FAST5_TO_POD5                     } from '../modules/pod5.nf'
include { BASECALL; DEMUX; DORADO_SUMMARY           } from '../modules/dorado.nf'
include { CONVERT_TO_FASTQ; MERGE_BAMS_FOR_SUMMARY  } from '../modules/samtools.nf'
include { PYCOQC                                    } from '../modules/pycoqc.nf'

def validateSingleFormat(listOfFormats){
    if (listOfFormats.size() != 1) {
        log.error("Multiple signal filetypes ${listOfFormats} found in '${params.raw_read_dir}'. Please separate filetypes into distinct directories and process indepedently.")
    }
}

workflow BASECALLING {  
    take:
    raw_read_signal_files

    main:
    //todo add reference aligned basecalling?

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
    | map{ barcode_kit, long_read_bam -> 
        def meta = [:]
        meta.barcode_kit = barcode_kit
        meta.barcode = "${ long_read_bam.simpleName.contains("barcode") ? long_read_bam.simpleName.split("barcode")[1] : long_read_bam.simpleName }" //i.e. when simpleName = unclassified
        tuple(meta, long_read_bam)
    }
    | set{ barcode_bam_ch }
    
    if (params.read_format == "fastq") {
        CONVERT_TO_FASTQ(barcode_bam_ch)
        | set { long_reads_ch }
    } else {
        barcode_bam_ch.set { long_reads_ch }
    }

    if (params.qc_reads) {
        LONG_READ_QC(barcode_bam_ch)
        | set { sequencing_summary }
    } else {
        sequencing_summary = Channel.of('')
    }

    emit:
    long_reads_ch
    sequencing_summary
}

workflow LONG_READ_QC {
    take:
    barcode_bam_ch
    
    main:
    
    barcode_bam_ch.map { meta, long_bam -> long_bam }
        | filter { long_bam -> long_bam.name != "unclassified.bam"}
        | collect
        | MERGE_BAMS_FOR_SUMMARY
        | DORADO_SUMMARY
        | PYCOQC
    
    emit:
    DORADO_SUMMARY.out.summary_channel
}