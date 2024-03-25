import org.apache.commons.io.FilenameUtils

include { CONVERT_FAST5_TO_POD5                                                 } from '../modules/pod5.nf'
include { MODEL_DOWNLOAD; BASECALL; DEMUX; DORADO_SUMMARY                       } from '../modules/dorado.nf'
include { CONVERT_TO_FASTQ; MERGE_BAMS_FOR_SUMMARY; REMOVE_DUPLICATES_FROM_BAMS } from '../modules/samtools.nf'
include { PYCOQC                                                                } from '../modules/pycoqc.nf'

def validateSingleFormat(listOfFormats){
    if (listOfFormats.size() != 1) {
        log.error("Multiple signal filetypes ${listOfFormats} found in '${params.raw_read_dir}'. Please separate filetypes into distinct directories and process indepedently.")
    }
}

def mark_read_duplicates_in_summary(sequencing_summary, outputFilePath){
    def filePath = sequencing_summary.toString()

    //make a directory in work if it doesn't exist
    def baseDir = new File(outputFilePath)

    if(!baseDir.exists()) {
        baseDir.mkdir()
    }

    // Read the TSV file
    def tsvFile = new File(filePath)

    //make a set to store what we have seen
    def column2Set = new HashSet<String>()
    def duplicates=[]

    // Iterate over each line in the file
    tsvFile.eachLine { line ->
        // Split the line by tabs
        def columns = line.split('\t')
    
        // Ensure the line has at least two columns
        if (columns.size() >= 2) {
            def value = columns[1]
            if (column2Set.contains(value)) {
                duplicates << value
            } else {
                column2Set.add(value)
            }
        }
    }

    def finalPath = "${baseDir}/duplicates_${workflow.start}.txt"

    def outputFile = new File(finalPath)

    duplicates.each { value ->
        outputFile.append(value + '\n')
    }

    return finalPath
}

workflow BASECALLING {  
    take:
    raw_read_signal_files
    additional_metadata

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
    | MODEL_DOWNLOAD
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

    LONG_READ_QC(barcode_bam_ch)
    | set { sequencing_summary }

    if (params.barcode_kit_name.size() >= 2) {
        sequencing_summary.map{ mark_read_duplicates_in_summary(it, "${workflow.workDir}/summary_duplicates/") }
        | set { dupliate_read_list }

        barcode_bam_ch.combine(dupliate_read_list)
        | REMOVE_DUPLICATES_FROM_BAMS
        | set { bam_ch }

    } else {
        barcode_bam_ch.set { bam_ch }
    }
    
    bam_ch
    | map { meta, reads -> ["${meta.barcode_kit}_${meta.barcode}", meta, reads]}
    | set { bam_by_barcode_ch }

    additional_metadata
    | join(bam_by_barcode_ch)
    | map { barcodekit_barcode, meta1, meta2, reads -> [meta1 + meta2, reads] }
    | set { bam_with_metadata_ch }

    if (params.read_format == "fastq") {
        CONVERT_TO_FASTQ(bam_with_metadata_ch)
        | set { long_reads_ch }

    } else {
        bam_ch.set { long_reads_ch }
    }

    emit:
    long_reads_ch
    sequencing_summary
    //model_ch = MODEL_DOWNLOAD.out.model_ch.map{ pod5, model -> model} //currently not useful as we use 9.4.1 flow cells but later this could be useful
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