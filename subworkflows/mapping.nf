include { SAMTOOLS_SORT; CONVERT_TO_BAM; SAMTOOLS_INDEX_BAM } from '../modules/samtools.nf'
include { MINIMAP2_INDEX; MINIMAP2_ALIGN } from '../modules/minimap2.nf'

workflow MAPPING {
    take:
    reference
    fastqs

    emit:
    sorted_reads_bam

    main:
    MINIMAP2_INDEX(
        reference
    )
    | set{ ref_mm2_index }
    
    fastqs.combine(ref_mm2_index)
        .set { minimap2_input }

    MINIMAP2_ALIGN(
        minimap2_input
    )
    | CONVERT_TO_BAM
    | SAMTOOLS_SORT
    | SAMTOOLS_INDEX_BAM

    SAMTOOLS_INDEX_BAM.out.bam_index
        .set{ sorted_reads_bam }
} 