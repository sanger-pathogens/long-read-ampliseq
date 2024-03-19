include { CUT_PRIMERS } from '../modules/cutadapt.nf'

workflow PROCESS_FILTER_READS {
    take:
    fastqs

    main:
    Channel.fromPath(params.primers)
    | set { primers }

    fastqs
    | combine(primers)
    | set { cutadapt_input }

    cutadapt_input.view()
    
    CUT_PRIMERS(cutadapt_input)
    CUT_PRIMERS.out.trimmed_reads
    | set { trimmed_reads }

    emit:
    trimmed_reads
}
