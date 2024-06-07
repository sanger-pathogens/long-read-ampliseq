include {
    GET_READLENGTH_DISTRIBUTION;
    SAMTOOLS_DEPTH;
    SAMTOOLS_STATS;
    ON_AND_OFF_TARGET_STATS
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
    mapped_reads_bam

    main:
    Channel.value("post_map_qc")
    | set { qc_stage }

    SAMTOOLS_STATS(
        mapped_reads_bam,
        qc_stage
    )
    SAMTOOLS_STATS.out.stats_ch.set { ch_samtools_stats }

    emit:
    ch_samtools_stats
}

workflow POST_FILTER_QC {
    take:
    on_target_reads_bam
    off_target_reads_bam
    target_regions_bed

    main:

    Channel.value("post_filter_qc")
    | set { qc_stage }

    on_target_reads_bam.map { meta, bam, bai -> [meta, bam]}
    | set { samtools_stats_input }
    SAMTOOLS_STATS(
        samtools_stats_input,
        qc_stage
    )
    SAMTOOLS_STATS.out.stats_ch.set { ch_samtools_stats }

    GET_READLENGTH_DISTRIBUTION(
        on_target_reads_bam,
        qc_stage
    )

    COVERAGE_QC(
        on_target_reads_bam,
        target_regions_bed,
        qc_stage
    )

    ON_AND_OFF_TARGET_STATS(
        on_target_reads_bam.join(off_target_reads_bam)
    )

    ON_AND_OFF_TARGET_STATS.out
        .collectFile(name: "${params.outdir}/qc/post_filter_qc/on_and_off_target_stats.csv", keepHeader: true, skip: 1) { it[1] }
        .set { ch_on_and_off_target_stats }

    emit:
    ch_samtools_stats
    ch_on_and_off_target_stats
}

workflow COVERAGE_QC {
    take:
    sorted_reads_bam
    target_regions_bed
    qc_stage

    main:
    BEDTOOLS_GENOMECOV(
        sorted_reads_bam,
        qc_stage
    )

    sorted_reads_bam.combine(target_regions_bed)
        .set { bedtools_coverage_input }
    
    BEDTOOLS_COVERAGE(
        bedtools_coverage_input,
        qc_stage
    )

    SAMTOOLS_DEPTH(
        sorted_reads_bam,
        qc_stage
    )

    SAMTOOLS_DEPTH.out.samtools_coverage
        .combine(target_regions_bed)
        .set { coverage_over_defined_regions_input }

    PYTHON_COVERAGE_OVER_DEFINED_REGIONS(
        coverage_over_defined_regions_input,
        qc_stage
    )

    PYTHON_COVERAGE_OVER_DEFINED_REGIONS.out.coverage_summary
        .collect() { it[1] }
        .set { coverage_summaries }

    PYTHON_PLOT_COVERAGE(
        coverage_summaries,
        qc_stage
    )
}