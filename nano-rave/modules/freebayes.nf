process FREEBAYES_VARIANT_CALLING {
    container "docker.io/gfanz/freebayes@sha256:d32bbce0216754bfc7e01ad6af18e74df3950fb900de69253107dc7bcf4e1351"
    input:
        tuple val(ref_id), path(reference)
        tuple path(sorted_bam_file), path(sorted_bam_index)
    output:
        path("*.vcf"), emit: vcf_ch
    script:
        """
        filename=\$(basename ${sorted_bam_file} | awk -F "." '{ print \$1}')
        
        freebayes -f ${reference} ${sorted_bam_file} > \${filename}.vcf
        
        # sort vcf by index to stop tabix crying
        cat \${filename}.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > \${filename}_sorted.vcf
        mv \${filename}_sorted.vcf \${filename}.vcf
        """
}