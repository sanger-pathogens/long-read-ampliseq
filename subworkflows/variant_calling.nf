include { CLAIR3_CALL } from '../modules/clair3.nf'


workflow CALL_VARIANTS {
    take:
    sorted_reads_bam_with_bed
    reference_index_ch

    main:
    sorted_reads_bam_with_bed.combine(reference_index_ch)
    | CLAIR3_CALL

    emit:
    CLAIR3_CALL.out.clair3_out
}