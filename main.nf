#!/usr/bin/env nextflow

//
// SUBWORKFLOWS
//

include { BASECALLING } from './subworkflows/basecalling.nf'
//include { NANO_RAVE   } from './nano-rave/subworkflow/nano-rave.nf'
include { PRE_MAP_QC } from './subworkflows/pre_map_qc.nf'
include { MAPPING } from './subworkflows/mapping.nf'

def logo = NextflowTool.logo(workflow, params.monochrome_logs)

log.info logo


def printHelp() {
    NextflowTool.help_message("${workflow.ProjectDir}/schema.json",
    params.monochrome_logs, log)
}

workflow {
    if (params.help) {
        printHelp()
        exit 0
    }
    if (params.basecall) {
        raw_reads = Channel.fromPath("${params.raw_read_dir}/*.{fast5,pod5}", checkIfExists: true)
        BASECALLING(raw_reads)
    }
    
    Channel.fromPath(params.reference)
        .set{ reference }

    PRE_MAP_QC(
        BASECALLING.out.long_reads_ch
    )

    MAPPING(
        reference,
        BASECALLING.out.long_reads_ch
    )
}

workflow.onComplete {
    NextflowTool.summary(workflow, params, log)
}