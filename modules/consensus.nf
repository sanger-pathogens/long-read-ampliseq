process CURATE_CONSENSUS {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/curated_consensus", pattern: "*.fasta", mode: 'copy', overwrite: true

    conda 'conda-forge::python=3.10.2'
    container 'quay.io/sangerpathogens/pysam:0.0.2'

    input:
    tuple val(meta), file(vcf_final), path(reference), path(ref_index)

    output:
    tuple val(meta), path("${meta.ID}.fasta"),  emit: multi_fasta
    tuple val(meta), path("${meta.ID}_multi_locus.fasta"), emit: full_consensus
    path("*.fa"), emit: per_loci

    script:
    """
    vcf_to_fasta.py -r ${reference} -v ${vcf_final} --fasta_id ${meta.ID} -b ${params.target_regions_bed} -rr --multi_locus --multifasta -o ${meta.ID}
    vcf_to_fasta.py -r ${reference} -v ${vcf_final} --fasta_id ${meta.ID} -b ${params.target_regions_bed} --singlefasta -o ${meta.ID}
    """

}