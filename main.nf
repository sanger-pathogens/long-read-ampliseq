#!/usr/bin/env nextflow

include { INDEX_REF } from './modules/samtools.nf'
include { MULTIQC } from './modules/multiqc.nf'

//
// SUBWORKFLOWS
//

include { BASECALLING } from './subworkflows/basecalling.nf'
include {
    PRE_MAP_QC as PRE_MAP_QC_PRE_TRIM;
    PRE_MAP_QC as PRE_MAP_QC_POST_TRIM;
} from './subworkflows/pre_map_qc.nf'
include { PROCESS_FILTER_READS } from './subworkflows/process_filter_reads.nf'
include { MAPPING } from './subworkflows/mapping.nf'
include { FILTER_BAM } from './subworkflows/post_map_filtering.nf'
include {
    POST_MAP_QC;
    POST_FILTER_QC;
} from './subworkflows/post_map_qc.nf'
include { CALL_VARIANTS } from './subworkflows/variant_calling.nf'

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
        .ifEmpty {exit 1, "${params.additional_metadata} appears to be an empty file!"}
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

    INDEX_REF(reference)
    | set { reference_index_ch }

    PRE_MAP_QC_PRE_TRIM(
        BASECALLING.out.long_reads_ch
    )

    BASECALLING.out.long_reads_ch.filter{ meta, bam -> meta.barcode != "unclassified"}
    | set { remove_unclassified_for_mapping }

    PROCESS_FILTER_READS(
        remove_unclassified_for_mapping
    )

    PRE_MAP_QC_POST_TRIM(
        PROCESS_FILTER_READS.out.trimmed_reads
    )

    MAPPING(
        reference,
        PROCESS_FILTER_READS.out.trimmed_reads
    )

    FILTER_BAM(
        MAPPING.out.mapped_reads_bam,
        target_regions_bed
    )

    POST_FILTER_QC(
        FILTER_BAM.out.on_target_reads_bam,
        FILTER_BAM.out.off_target_reads_bam,
        target_regions_bed
    )

    //for now until we decide filter logic we can just do some basic variant calling (large unfiltered data though)
    CALL_VARIANTS(
        FILTER_BAM.out.on_target_reads_bam.combine(target_regions_bed),
        reference_index_ch
    )

    MULTIQC(
        BASECALLING.out.pycoqc_json.ifEmpty([]),
        PRE_MAP_QC_PRE_TRIM.out.ch_fastqc_raw_zip.collect{it[1]}.ifEmpty([]),
        PRE_MAP_QC_POST_TRIM.out.ch_fastqc_raw_zip.collect{it[1]}.ifEmpty([]),
        MAPPING.out.ch_samtools_stats.collect{it[1,2]}.ifEmpty([]),
        POST_FILTER_QC.out.ch_samtools_stats.collect{it[1,2]}.ifEmpty([])
    )
}

workflow.onComplete {
    NextflowTool.summary(workflow, params, log)
}
