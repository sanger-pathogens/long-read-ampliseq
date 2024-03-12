process GET_CHROM_SIZES_AND_INDEX {
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    input:
        tuple val(ref_id), path(reference)
    output:
        tuple val(ref_id), path("*.fai"), emit: fasta_index_ch
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        samtools faidx ${reference}
        """
}