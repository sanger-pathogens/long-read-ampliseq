process MINIMAP2_INDEX {
    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"
    input:
        tuple val(ref_id), path(reference)
    output:
        path("*.mmi"), emit: mm2_index_ch
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        minimap2 -ax map-ont -t 2 -d ${reference}.mmi ${reference}
        """
}

process MINIMAP2_ALIGN {
    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"
    input:
        tuple val(ref_id), path(reference)
        path(mm2_index)
        each path(fastq)
    output:
        path("*.sam"), emit: sam_file
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        ref_id=\$(echo $reference | awk -F "seq_" '{ print \$NF }' | sed 's|.fasta||g')
        fname=\$(basename $fastq | awk -F "." '{ print \$1}')
        minimap2 -ax map-ont --MD -t 2 $mm2_index $fastq > \${fname}_\${ref_id}.sam
        """
}