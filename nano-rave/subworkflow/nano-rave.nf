#!/usr/bin/env nextflow
// Copyright (C) 2022,2023 Genome Research Ltd.

include { BEDTOOLS_GENOMECOV                                                           } from '../modules/bedtools.nf'
include { SORT_FASTQS; NORMALISE_FASTAS                                                } from '../modules/file_utility.nf'
include { MEDAKA_VARIANT_CALLING; MEDAKA_HAPLOID_VARIANT_CALLING                       } from '../modules/medaka.nf'
include { NANOPLOT_QC                                                                  } from '../modules/nanoplot.nf'
include { GET_CHROM_SIZES_AND_INDEX; SAMTOOLS_VIEW_SAM_TO_BAM; SAMTOOLS_SORT_AND_INDEX } from '../modules/samtools.nf'
include { CLAIR3_VARIANT_CALLING                                                       } from '../modules/clair3.nf'
include { FREEBAYES_VARIANT_CALLING                                                    } from '../modules/freebayes.nf'
include { MINIMAP2_INDEX; MINIMAP2_ALIGN                                               } from '../modules/minimap2.nf'
include { PYCOQC                                                                       } from '../modules/pycoqc.nf'
include { BGZIP_AND_INDEX_VCF                                                          } from '../modules/tabix.nf'


def validate_path_param(
    param_option, 
    param, 
    type="file", 
    mandatory=true) 
{
    valid_types=["file", "directory"]
    if (!valid_types.any { it == type }) {
            log.error("Invalid type '${type}'. Possibilities are ${valid_types}.")
            return 1
    }
    param_name = (param_option - "--").replaceAll("_", " ")
    if (param) {
        def file_param = file(param)
        if (!file_param.exists()) {
            log.error("The given ${param_name} '${param}' does not exist.")
            return 1
        } else if (
              (type == "file" && !file_param.isFile())
              ||
              (type == "directory" && !file_param.isDirectory())
          ) {
            log.error("The given ${param_name} '${param}' is not a ${type}.")
            return 1
        }
    } else if (mandatory) {
        log.error("No ${param_name} specified. Please specify one using the ${param_option} option.")
        return 1
    }
    return 0
}

def validate_choice_param(param_option, param, choices) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (param) {
        if (!choices.any { it.contains(param.toString()) }) {
            log.error("Please specify the ${param_name} using the ${param_option} option. Possibilities are ${choices}.")
            return 1
        }
    } else {
        log.error("Please specify the ${param_name} using the ${param_option} option")
        return 1
    }
    return 0
}

def validate_number_param(param_option, param) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (param != null) /* Explicit comparison with null, because 0 is an acceptable value */ {
        if (!(param instanceof Number)) {
            log.error("The ${param_name} specified with the ${param_option} option must be a valid number")
            return 1
        }
    } else {
        log.error("Please specify the ${param_name} using the ${param_option} option")
        return 1
    }
    return 0
}

def validate_min_barcode_dir_size(param_option, param) {
    if (validate_number_param(param_option, param) == 1) {
        return 1
    }
    param_name = (param_option - "--").replaceAll("_", " ")
    if (!(param > 0)) {
        log.error("The ${param_name} specified with the ${param_option} option must have a positive value")
        return 1
    }
    return 0
}

def validate_results_dir(results_dir) {
    results_dir = file(results_dir)
    if (results_dir.exists() && !results_dir.isDirectory()) {
        log.error("The given results_dir '${results_dir}' is not a directory.")
        return 1
    }
    return 0
}

def validate_clair3_args(clair3_args) {
    def invalid_options = ['--bam_fn', '--ref_fn', '--threads', '--platform', '--output']
    def required_options = ['--model_path']
    def errors = 0
    if (params.variant_caller != "clair3" && clair3_args) {
        log.error("Clair3 arguments were provided but clair3 was not set as the --variant_caller!")
        errors += 1
    } else if (params.variant_caller == "clair3" && clair3_args) {
        def invalid_options_found = invalid_options.findAll { clair3_args.contains(it) }
        if (invalid_options_found) {
            log.error("The following clair3 options were provided in --clair3_args, but are reserved for use by this pipeline: ${invalid_options_found}")
            errors += 1
        }
        def required_options_found = required_options.findAll { clair3_args.contains(it) }
        def required_options_not_found = required_options.minus(required_options_found)
        if (required_options_not_found) {
            log.error("The following clair3 options were not provided in --clair3_args, but are required for use by this pipeline: ${required_options_not_found}")
            errors += 1
        }
    }
    return errors
}

def validate_parameters() {
    def errors = 0

    errors += validate_path_param("--reference_manifest", params.reference_manifest)
    errors += validate_path_param("--sequencing_manifest", params.sequencing_manifest)
    errors += validate_choice_param("--variant_caller", params.variant_caller, ["medaka", "medaka_haploid", "freebayes", "clair3"])
    errors += validate_min_barcode_dir_size("--min_barcode_dir_size", params.min_barcode_dir_size)
    errors += validate_results_dir(params.results_dir)
    errors += validate_clair3_args(params.clair3_args)

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}

workflow NANOSEQ {
    take:
        fastq_ch

    main:
        ref_manifest_ch = Channel.fromPath(params.reference_manifest)
        ref_path_ch = ref_manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.reference_id, file(row.reference_path)) }

        NANOPLOT_QC(fastq_ch)

        NORMALISE_FASTAS(ref_path_ch)

        GET_CHROM_SIZES_AND_INDEX(
            NORMALISE_FASTAS.out.normalised_ref_ch
        )

        MINIMAP2_INDEX(
            GET_CHROM_SIZES_AND_INDEX.out.ref_ch
        )

        MINIMAP2_ALIGN(
            MINIMAP2_INDEX.out.ref_ch,
            MINIMAP2_INDEX.out.mm2_index_ch,
            fastq_ch
        )

        SAMTOOLS_VIEW_SAM_TO_BAM(
            MINIMAP2_ALIGN.out.ref_ch,
            MINIMAP2_ALIGN.out.sam_file
        )

        SAMTOOLS_SORT_AND_INDEX(
            SAMTOOLS_VIEW_SAM_TO_BAM.out.ref_ch,
            SAMTOOLS_VIEW_SAM_TO_BAM.out.bam_file,
        )

        BEDTOOLS_GENOMECOV(
            SAMTOOLS_SORT_AND_INDEX.out.sorted_bam
        )

        if (params.variant_caller == "medaka_haploid") {
            MEDAKA_HAPLOID_VARIANT_CALLING(
                NORMALISE_FASTAS.out.normalised_ref_ch,
                fastq_ch
            )
            BGZIP_AND_INDEX_VCF(MEDAKA_HAPLOID_VARIANT_CALLING.out.vcf_ch)
        }

        if (params.variant_caller == "medaka") {
            MEDAKA_VARIANT_CALLING(
                SAMTOOLS_SORT_AND_INDEX.out.ref_ch,
                SAMTOOLS_SORT_AND_INDEX.out.sorted_bam
            )
            BGZIP_AND_INDEX_VCF(MEDAKA_VARIANT_CALLING.out.vcf_ch)
        }

        if (params.variant_caller == "freebayes") {
            FREEBAYES_VARIANT_CALLING(
                SAMTOOLS_SORT_AND_INDEX.out.ref_ch,
                SAMTOOLS_SORT_AND_INDEX.out.sorted_bam
            )
            BGZIP_AND_INDEX_VCF(FREEBAYES_VARIANT_CALLING.out.vcf_ch)
        }

        if (params.variant_caller == "clair3") {
            SAMTOOLS_SORT_AND_INDEX.out.ref_ch
                .combine(GET_CHROM_SIZES_AND_INDEX.out.fasta_index_ch, by: 0)
                .set { clair3_ref_input }
            CLAIR3_VARIANT_CALLING(
                clair3_ref_input,
                SAMTOOLS_SORT_AND_INDEX.out.sorted_bam
            )
            BGZIP_AND_INDEX_VCF(CLAIR3_VARIANT_CALLING.out.vcf_ch)
        }
}

workflow NANO_RAVE {
    take:
    sequence_ch

    /*
    validate_parameters()

    sequencing_manifest = Channel.fromPath(sequencing_manifest)
    sequence_ch = sequencing_manifest.splitCsv(header: true, sep: ',')
                    .map{ row -> tuple(file(row.sequencing_dir), file(row.sequence_summary_file)) }

    */

    SORT_FASTQS(sequence_ch)

    /*
    excluded as it is run in the barcoding process

    PYCOQC(sequence_ch)
    */

    NANOSEQ(SORT_FASTQS.out.full_fastq_files.flatten())
}