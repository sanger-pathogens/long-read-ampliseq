include { FASTQC } from '../modules/fastqc.nf'

workflow PRE_MAP_QC {
    take:
    fastqs

    main:
    FASTQC(
        fastqs
    )
    FASTQC.out.zip.set { ch_fastqc_raw_zip }
    FASTQC.out.html.set { ch_fastqc_raw_html }

    emit:
    ch_fastqc_raw_zip
    ch_fastqc_raw_html
}
