include {
    REMOVE_OFF_TARGET_READS;
    SAMTOOLS_INDEX_BAM as INDEX_ON_TARGET_BAM;
    SAMTOOLS_INDEX_BAM as INDEX_OFF_TARGET_BAM;
} from '../modules/samtools.nf'

workflow FILTER_BAM {
    take:
    mapped_reads_bam
    target_regions_bed

    main:
    mapped_reads_bam
    | combine(target_regions_bed)
    | set { remove_off_targets_input }

    REMOVE_OFF_TARGET_READS(remove_off_targets_input)
    INDEX_ON_TARGET_BAM(REMOVE_OFF_TARGET_READS.out.on_target_reads_bam)
    INDEX_OFF_TARGET_BAM(REMOVE_OFF_TARGET_READS.out.off_target_reads_bam)

    INDEX_ON_TARGET_BAM.out.bam_index
    | set { on_target_reads_bam }
    INDEX_OFF_TARGET_BAM.out.bam_index
    | set { off_target_reads_bam }

    emit:
    on_target_reads_bam
    off_target_reads_bam
}
