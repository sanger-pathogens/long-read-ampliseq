process BGZIP_AND_INDEX_VCF {
    container "quay.io/biocontainers/tabix:1.11--hdfd78af_0"
    publishDir "${params.results_dir}/variant_calling", mode: 'copy', overwrite: true
    input:
        path(vcf_file)

    output:
        path("*.vcf.gz"), emit: bgzip_vcf_file_ch
        path("*.vcf.gz.tbi"), emit: vcf_index_ch

    script:
        """
        filename=\$(basename $vcf_file| awk -F "." '{ print \$1}')
        bgzip -c $vcf_file > \${filename}.vcf.gz
        tabix \${filename}.vcf.gz
        """
}