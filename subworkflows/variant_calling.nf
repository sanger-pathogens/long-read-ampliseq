include { CLAIR3_CALL } from '../modules/clair3.nf'
include { GUNZIP } from '../modules/helper_process.nf'
include { CURATE_CONSENSUS } from '../modules/consensus.nf'
include { MERGE_GVCF; BCFTOOLS_QUERY } from '../modules/bcftools.nf'
include { CONSTRUCT_PHYLO } from '../assorted-sub-workflows/tree_build/tree_build.nf'

workflow CALL_VARIANTS {
    take:
    filtered_reads_bam_with_bed
    reference_index_ch

    main:
    filtered_reads_bam_with_bed
    | combine(reference_index_ch)
    | combine(Channel.fromPath(params.clair3_model))
    | CLAIR3_CALL

    GUNZIP( CLAIR3_CALL.out.vcf_out )
    | combine(Channel.fromPath(params.target_regions_bed))
    | CURATE_CONSENSUS

    CLAIR3_CALL.out.clair3_out.map{ metadata, path -> path }
    | collect
    | MERGE_GVCF
    | BCFTOOLS_QUERY

    CURATE_CONSENSUS.out.full_consensus.collectFile { meta, file -> [ "merged.fasta", file ] }
    | CONSTRUCT_PHYLO

    //CURATE_CONSENSUS.out.per_loci.flatten().map{ loci_fasta -> 
    //    
    //}

    emit:
    CLAIR3_CALL.out.clair3_out
}
