include { GET_READLENGTH_DISTRIBUTION } from '../modules/samtools.nf'
include { BEDTOOLS_GENOMECOV; BEDTOOLS_COVERAGE } from '../modules/bedtools.nf'

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
}
