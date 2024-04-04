include { REMOVE_OFF_TARGET_READS; SAMTOOLS_INDEX_BAM } from '../modules/samtools.nf'

workflow FILTER_BAM {
    take:
    sorted_reads_bam
    target_regions_bed

    main:

    sorted_reads_bam
    | combine(target_regions_bed)
    | set { remove_off_targets_input }

    
    REMOVE_OFF_TARGET_READS(remove_off_targets_input)

    SAMTOOLS_INDEX_BAM(REMOVE_OFF_TARGET_READS.out.on_target_reads_bam)
    SAMTOOLS_INDEX_BAM.out.bam_index
    | set { on_target_reads_bam }

    REMOVE_OFF_TARGET_READS.out.off_target_reads_bam
    | set { off_target_reads_bam }

    emit:
    on_target_reads_bam
    off_target_reads_bam
}
