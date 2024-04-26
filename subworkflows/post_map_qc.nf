include {
    GET_READLENGTH_DISTRIBUTION;
    SAMTOOLS_DEPTH;
    SAMTOOLS_STATS;
    COUNT_ON_AND_OFF_TARGET_READS
} from '../modules/samtools.nf'
include {
    BEDTOOLS_GENOMECOV;
    BEDTOOLS_COVERAGE
} from '../modules/bedtools.nf'
include {
    PYTHON_COVERAGE_OVER_DEFINED_REGIONS;
    PYTHON_PLOT_COVERAGE
} from '../modules/custom.nf'

workflow POST_MAP_QC {
    take:
    sorted_reads_bam
    on_target_reads_bam
    off_target_reads_bam
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

    SAMTOOLS_DEPTH.out.samtools_coverage
        .combine(target_regions_bed)
        .set { coverage_over_defined_regions_input }

    PYTHON_COVERAGE_OVER_DEFINED_REGIONS(
        coverage_over_defined_regions_input
    )

    SAMTOOLS_STATS(
        sorted_reads_bam
    )

    COUNT_ON_AND_OFF_TARGET_READS(
        on_target_reads_bam.join(off_target_reads_bam)
    )

    PYTHON_COVERAGE_OVER_DEFINED_REGIONS.out.coverage_summary
        .collect()
        .set { coverage_summaries }

    PYTHON_PLOT_COVERAGE(coverage_summaries)
}
