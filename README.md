# Long read AmpliSeq

## Installation
1. [Install Nextflow](https://www.nextflow.io/docs/latest/install.html)
2. [Install Docker](https://docs.docker.com/engine/install/)
3. Clone the repository

    Please note you will also need access to the [assorted-sub-workflows](https://github.com/sanger-pathogens/assorted-sub-workflows/tree/b065b17b0ee663483fa14c09fc9b1dede9afa8ba) and [nextflowtool](https://github.com/sanger-pathogens/nextflowtool/tree/74b25a9346d243db662caccb777296061400b65a) submodule repos

    With SSH (will need an SSH key- instructions [here](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)):
    ```
    git clone --recurse-submodules git@github.com:sanger-pathogens/long-read-ampliseq.git
    ```

    With HTTPS:
    
    ```
    git clone --recurse-submodules https://github.com/sanger-pathogens/long-read-ampliseq.git
    ```
    You may need to enter a personal access token for authentication in place of a password- please see this [guide](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens). 

## Usage
```
nextflow run long-read-ampliseq/main.nf \
--raw_read_dir <directory containing FAST5/POD5 files> \
--reference <reference fasta> \
--primers <fasta containing primers> \
--target_regions_bed <BED file containing target regions> \
--additional_metadata <CSV mapping sample IDs to barcodes> \
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
- --remove_recombination = false
- --raxml_base_model = 'GTR+G4'
- --raxml_threads = 2


## Support
Please contact PaM Informatics for support through our [helpdesk portal](https://jira.sanger.ac.uk/servicedesk/customer/portal/16) or for external users please reach out by email: pam-informatics@sanger.ac.uk
