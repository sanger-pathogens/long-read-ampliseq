process CLAIR3_VARIANT_CALLING {
    container "quay.io/biocontainers/clair3:1.0.5--py39hf5e1c6e_0"  // Includes models
    input:
        tuple val(ref_id), path(reference), path(reference_index)
        tuple path(sorted_bam_file), path(sorted_bam_index)
    output:
        path("*.vcf"), emit: vcf_ch
    script:
        """
        filename=\$(basename ${sorted_bam_file} | awk -F "." '{ print \$1}')
        
        run_clair3.sh \
            --bam_fn=${sorted_bam_file} \
            --ref_fn=${reference} \
            --threads=${task.cpus} \
            --platform="ont" \
            --model_path="/opt/models/${params.clair3_model}"
            --output=. \
            ${params.clair3_args}

        if [[ ! -s "merge_output.vcf.gz" ]]; then
            cp pileup.vcf.gz merge_output.vcf.gz
        fi
        # sort vcf by index to stop tabix crying
        zcat merge_output.vcf.gz | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > \${filename}.vcf
        """
}