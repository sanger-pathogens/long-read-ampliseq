include {
    SAMTOOLS_SORT;
    CONVERT_TO_BAM;
    SAMTOOLS_INDEX_BAM;
} from '../modules/samtools.nf'
include {
    MINIMAP2_INDEX;
    MINIMAP2_ALIGN;
} from '../modules/minimap2.nf'

include { POST_MAP_QC } from './post_map_qc.nf'

workflow MAPPING {
    take:
    reference
    fastqs

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

    POST_MAP_QC(
        MINIMAP2_ALIGN.out.mapped_reads_sam,
    )
    POST_MAP_QC.out.ch_samtools_stats
    | set { ch_samtools_stats }

    CONVERT_TO_BAM(
        MINIMAP2_ALIGN.out.mapped_reads_sam
    )
    | SAMTOOLS_SORT
    | SAMTOOLS_INDEX_BAM
    | set { mapped_reads_bam }

    emit:
    mapped_reads_bam
    ch_samtools_stats
} 