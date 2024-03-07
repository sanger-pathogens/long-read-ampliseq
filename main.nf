#!/usr/bin/env nextflow

//
// SUBWORKFLOWS
//

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
}