include { GET_READLENGTH_DISTRIBUTION; SAMTOOLS_DEPTH } from '../modules/samtools.nf'
include { BEDTOOLS_GENOMECOV; BEDTOOLS_COVERAGE } from '../modules/bedtools.nf'
include { PYTHON_COVERAGE_OVER_DEFINED_REGIONS } from '../modules/custom.nf'

workflow POST_MAP_QC {
    take:
    sorted_reads_bam
    target_regions_bed

    main:
    GET_READLENGTH_DISTRIBUTION(
        sorted_reads_bam
    )

    BEDTOOLS_GENOMECOV(
        sorted_reads_bam
    )

    sorted_reads_bam.combine(target_regions_bed)
        .set { bedtools_coverage_input }
    

    BEDTOOLS_COVERAGE(
        bedtools_coverage_input
    )

    SAMTOOLS_DEPTH(
        sorted_reads_bam
    )

    sorted_reads_bam
        .join(SAMTOOLS_DEPTH.out.samtools_coverage)
        .combine(target_regions_bed)
        .set { coverage_over_defined_regions_input }

    PYTHON_COVERAGE_OVER_DEFINED_REGIONS(
        coverage_over_defined_regions_input
    )
}