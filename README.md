# Long read AmpliSeq

## Installation
1. [Install Nextflow](https://www.nextflow.io/docs/latest/install.html)
2. Clone the repository

    With SSH:
    ```
    git clone --recurse-submodules git@gitlab.internal.sanger.ac.uk:sanger-pathogens/pipelines/long-read-ampliseq.git
    ```

    With HTTPS:
    ```
    git clone --recurse-submodules https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/long-read-ampliseq.git
    ```

## Usage
```
nextflow run long-read-ampliseq/main.nf \
--raw_read_dir <directory containing FAST5/POD5 files> \
--reference <reference fasta> \
--primers <fasta containing primers> \
--target_regions_bed <BED file containing target regions> \
--additional_metadata <CSV of additional metadata> \
-profile docker
```

##### Other parameters:

###### Basecalling
- --basecall = "true"
- --basecall_model = "dna_r9.4.1_e8_hac@v3.3"
- --trim_adapters = "all"
- --barcode_kit_name = ["EXP-NBD104", "EXP-NBD114"]
- --read_format = "fastq"

###### Saving output files
- --keep_sorted_bam = true
- --save_fastqs = true
- --save_trimmed = true
- --save_too_short = true
- --save_too_long = true

###### QC
- --qc_reads = true
- --min_qscore = 9
- --cutadapt_args = "-e 0.15 --no-indels --overlap 18"
- --lower_read_length_cutoff = 450
- --upper_read_length_cutoff = 800
- --coverage_reporting_thresholds = "1,2,8,10,25,30,40,50,100"
- --coverage_filtering_threshold = "25"

###### Variant calling
- --clair3_model = "r941_prom_hac_g360+g422"

## Support
Please contact PaM Informatics for support through our [helpdesk portal](https://jira.sanger.ac.uk/servicedesk/customer/portal/16) or for external users please reach out by email: pam-informatics@sanger.ac.uk
