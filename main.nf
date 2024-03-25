#!/usr/bin/env nextflow

//
// SUBWORKFLOWS
//

include { BASECALLING } from './subworkflows/basecalling.nf'
include { PRE_MAP_QC } from './subworkflows/pre_map_qc.nf'
include { MAPPING } from './subworkflows/mapping.nf'
include { POST_MAP_QC } from './subworkflows/post_map_qc.nf'
include { PROCESS_FILTER_READS } from './subworkflows/process_filter_reads.nf'

def logo = NextflowTool.logo(workflow, params.monochrome_logs)

log.info logo


def printHelp() {
    NextflowTool.help_message("${workflow.ProjectDir}/schema.json",
    params.monochrome_logs, log)
}

def combine_metadata_maps(barcode, meta1, meta2, reads) {
    [meta1 + meta2, reads]
}

workflow {
    if (params.help) {
        printHelp()
        exit 0
    }

    Channel.fromPath(params.reference)
        .set{ reference }

    Channel.fromPath(params.target_regions_bed)
        .set{ target_regions_bed }

    Channel.fromPath(params.additional_metadata)
        .ifEmpty {exit 1, "${params.additional_metedata} appears to be an empty file!"}
        .splitCsv(header:true, sep:',')
        .map { meta -> ["${meta.barcode_kit}_${meta.barcode}", meta] }
        .set { additional_metadata }

    if (params.basecall) {
        raw_reads = Channel.fromPath("${params.raw_read_dir}/*.{fast5,pod5}", checkIfExists: true)
        BASECALLING(
            raw_reads,
            additional_metadata
        )
    }

    PRE_MAP_QC(
        BASECALLING.out.long_reads_ch
    )

    MAPPING(
        reference,
        BASECALLING.out.long_reads_ch
    )

    POST_MAP_QC(
        MAPPING.out.sorted_reads_bam,
        target_regions_bed
    )

    PROCESS_FILTER_READS(
        BASECALLING.out.long_reads_ch
    )
}

workflow.onComplete {
    NextflowTool.summary(workflow, params, log)
}