process CLAIR3_CALL {
    tag "${meta.ID}"
    label "cpu_8"
    label "mem_32"
    label "time_12"

    container  'hkubal/clair3:v1.0.9'

    publishDir "${params.outdir}/variants/", mode: 'copy', overwrite: true, pattern: 'merge_output.gvcf.gz', saveAs: { filename -> "${meta.ID}_clair3.gvcf.gz" }
    publishDir "${params.outdir}/variants/logs/", mode: 'copy', overwrite: true, pattern: 'run_clair3.log', saveAs: { filename -> "${meta.ID}_clair3.log" }

    input:
    tuple val(meta), path(filtered_bam), path(bam_index), path(target_regions_bed), path(reference), path(reference_index)

    output:
    tuple val(meta), path("merge_output.gvcf.gz"), emit: clair3_out
    tuple val(meta), path("merge_output.vcf.gz"), path(reference), path(reference_index),  emit: vcf_out
    path("run_clair3.log")


    script:
    """
    run_clair3.sh \\
    --bam_fn=${filtered_bam} \\
    --ref_fn=${reference} \\
    --threads=${task.cpus} \\
    --platform="ont" \\
    --model_path="${params.clair3_model}" \\
    --output=. \\
    --sample_name=${meta.ID} \\
    --bed_fn=${target_regions_bed} \\
    --include_all_ctgs \\
    --haploid_precise \\
    --min_coverage=${params.clair3_min_coverage} \\
    --call_snp_only \\
    --print_ref_calls \\
    --gvcf \\
    --var_pct_full=1 \\
    --ref_pct_full=1 \\
    --var_pct_phasing=1 \\
    --no_phasing_for_fa
    """
}
