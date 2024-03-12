process MEDAKA_VARIANT_CALLING {
    container "quay.io/biocontainers/medaka:1.4.4--py38h130def0_0"
    input:
        tuple val(ref_id), path(reference)
        tuple path(sorted_bam_file), path(sorted_bam_index)
    output:
        path("*.vcf"), emit: vcf_ch
    script:
        """
        filename=\$(basename $sorted_bam_file | awk -F "." '{ print \$1}')
        
        medaka_variant -f ${reference} -i ${sorted_bam_file}

        mv medaka_variant/round_1.vcf \${filename}.vcf

        # sort vcf by index to stop tabix crying
        cat \${filename}.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > \${filename}_sorted.vcf
        mv \${filename}_sorted.vcf \${filename}.vcf
        """
}

process MEDAKA_HAPLOID_VARIANT_CALLING {
    container "quay.io/biocontainers/medaka:1.4.4--py38h130def0_0"
    input:
        tuple val(ref_id), path(reference)
        each path(fastq)
    output:
        path("*.vcf"), emit: vcf_ch
    script:
        """
        filename=\$(basename ${fastq} | awk -F '.' '{ print \$1}')_${ref_id}

        medaka_haploid_variant -r ${reference} -i ${fastq}

        mv medaka/medaka.annotated.vcf \${filename}.vcf

        # sort vcf by index to stop tabix crying
        cat \${filename}.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > \${filename}_sorted.vcf
        mv \${filename}_sorted.vcf \${filename}.vcf
        """
}