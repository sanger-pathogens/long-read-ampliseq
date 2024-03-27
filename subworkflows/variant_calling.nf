include { CLAIR3_CALL } from '../modules/clair3.nf'
include { GUNZIP } from '../modules/helper_process.nf'

include { CURATE_CONSENSUS } from '../assorted-sub-workflows/strain_mapper/modules/curate.nf'

workflow CALL_VARIANTS {
    take:
    sorted_reads_bam_with_bed
    reference_index_ch

    main:
    sorted_reads_bam_with_bed.combine(reference_index_ch)
    | CLAIR3_CALL

    GUNZIP( CLAIR3_CALL.out.vcf_out )
    | CURATE_CONSENSUS

    emit:
    CLAIR3_CALL.out.clair3_out
}