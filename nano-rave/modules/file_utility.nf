process SORT_FASTQS {
    /* Checks each barcode directory contains sufficient reads and concatenates fastq.gz files for further processing */
    input:
        tuple val(sequencing_dir), val(sequence_summary_file)
    output:
        path("*.fastq.gz"), emit: full_fastq_files
    script:
        """
        threshold=${params.min_barcode_dir_size}
        sample_name=\$(echo $sequencing_dir | awk -F "/" '{ print \$(NF-1) }')
        echo \$sample_name
        for dir in ${sequencing_dir}/fastq_pass/barcode*
        do
          barcode=\$(echo \${dir} | awk -F "/" '{ print \$NF }')
          disk_usage=\$(du --apparent-size -shm \$dir | awk '{ print \$1 }')
          if [ "\${disk_usage}" -ge "\${threshold}" ]
          then
            zcat \${dir}/*.fastq.gz > \${sample_name}_\${barcode}.fastq
          else
            echo "WARN: Skipping '\${barcode}' directory '\${dir}' as it contains <\${threshold}MB of fastq.gz files." >&2
          fi
        done
        gzip *.fastq
        """
}

process NORMALISE_FASTAS {
    container "quay.io/biocontainers/biopython:1.78"
    input:
        tuple val(ref_id), path(reference)
    output:
        tuple val(ref_id), path(reference), emit: normalised_ref_ch
    script:
        """
        #!/usr/bin/env python3
        from Bio import SeqIO
        import shutil

        records = SeqIO.parse("${reference}", "fasta")
        SeqIO.write(records, "${reference}.normalised", "fasta")
        shutil.move("${reference}.normalised", "${reference}")
        """
}