process CURATE_CONSENSUS {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/curated_consensus", mode: 'copy', overwrite: true

    conda 'conda-forge::python=3.10.2'
    container 'quay.io/sangerpathogens/pysam:0.0.2'

    input:
    tuple val(meta), file(vcf_final), path(reference), path(ref_index)

    output:
    tuple val(meta), path("*_multi_locus.fa"),  emit: multi_fasta
    tuple val(meta), path("*_full.fa"), emit: full_consensus

    script:
    """
    vcf_to_fasta.py -r ${reference} -v ${vcf_final} --fasta_id ${meta.ID} -b ${params.target_regions_bed} -rr true -o ${meta.ID}_full.fa
    vcf_to_fasta.py -r ${reference} -v ${vcf_final} --fasta_id ${meta.ID} -b ${params.target_regions_bed} -rr true --multifasta true -o ${meta.ID}_multi_locus.fa
    """
}