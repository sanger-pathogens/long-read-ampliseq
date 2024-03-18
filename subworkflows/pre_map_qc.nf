include { FASTQC } from '../modules/fastqc.nf'

workflow PRE_MAP_QC {
    take:
    fastqs

    main:
    FASTQC(
        fastqs
    )
}
