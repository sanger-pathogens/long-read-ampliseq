# Long-read Ampliseq

A Nextflow pipeline for basecalling, read mapping, QC, variant calling and analysis of nanopore multiplex amplicon data.

## Installation
1. [Install Nextflow](https://www.nextflow.io/docs/latest/install.html)

2. [Install Docker](https://docs.docker.com/engine/install/)

3. Download the appropriate Dorado installer from the [repo](https://github.com/nanoporetech/dorado#installation). The path to the executable will be ```<path to downloaded folder>/bin/dorado```

4. (Optional) Download the appropriate Dorado model from the [repo](https://github.com/nanoporetech/dorado/#available-basecalling-models)
    ```
    # Download all models
    dorado download --model all
    # Download particular model
    dorado download --model <model>
    ```
    If a pre-downloaded model path is not provided to the pipeline, the model specified by the `--basecall_model` parameter will be downloaded on the fly.

5. Download the appropriate Clair3 model from the [Rerio repo](https://github.com/nanoporetech/rerio?tab=readme-ov-file#clair3-models) (you will need Python3)
    
    First clone the repo:
    ```
    git clone https://github.com/nanoporetech/rerio
    ```
    This contains scripts to download the model(s) to ```clair3_models/<config>```
    ```
    #  Download all models
    python3 download_model.py --clair3
    #  Download particular model
    python3 download_model.py --clair3 clair3_models/<config>_model
    ```
    Each downloaded model can be found in the repo directory under ```clair3_models/<config>```

6. Clone the repository with required submodules
    
    ```
    git clone --recurse-submodules https://github.com/sanger-pathogens/long-read-ampliseq.git
    ```

## Usage
```
nextflow run long-read-ampliseq/main.nf \
--raw_read_dir <directory containing FAST5/POD5 files> \
--reference <reference fasta> \
--primers <fasta containing primers> \
--target_regions_bed <BED file containing target regions> \
--additional_metadata <CSV mapping sample IDs to barcodes> \
--dorado_local_path <absolute path to Dorado executable> \
--clair3_model <path to Clair3 model> \
-profile docker
```
The [examples](examples) folder contains some example files.

Instead of `-profile docker`, you can run the pipeline with `-profile laptop`. As well as enabling docker, the laptop profile allows the pipeline to be used offline by providing a local copy of a configuration file that is otherwise downloaded.

Should you need to run the pipeline offline, it is best to make use of pre-populated dependency caches. These can be created with any of the supported profiles (e.g. `-profile docker`) and involves running the pipeline once to completion. You will also need to provide a `--basecall_model_path` (see installation step 4)- the laptop profile includes a default local path for this, as well as `--clair3_model` and `--basecall_model_path`.

You can override the default paths using the command line parameters directly when invoking nextflow or supplying an additional config file in which these parameters are set, using the `-c my_custom.config` nextflow option.

### Other parameters:

#### Basecalling
- --basecall = "true"
- --basecall_model = "dna_r10.4.1_e8.2_400bps_hac@v4.3.0"
- --basecall_model_path = ""
- --trim_adapters = "all"
- --barcode_kit_name = ["SQK-NBD114-24"] (currently this can only be edited via the config file)
- --read_format = "fastq"

#### Saving output files
- --keep_sorted_bam = true
- --save_fastqs = true
- --save_trimmed = true
- --save_too_short = true
- --save_too_long = true

#### QC
- --qc_reads = true
- --min_qscore = 9
- --cutadapt_args = "-e 0.15 --no-indels --overlap 18"
- --lower_read_length_cutoff = 450
- --upper_read_length_cutoff = 800
- --coverage_reporting_thresholds = "1,2,8,10,25,30,40,50,100"
- --coverage_filtering_threshold = "25"
- --multiqc_config = ""

#### Variant calling
- --clair3_min_coverage = "5"
- --masking_quality = "15"

###### Consensus curation
- --min_ref_gt_qual = 1
- --min_alt_gt_qual = 1

###### Tree building
- --remove_recombination = false
- --raxml_base_model = 'GTR+G4'
- --raxml_threads = 2


## Running on Sanger farm

Load nextflow and singularity modules:
```bash
module load nextflow ISG/singularity
```

Clone the repository with required submodules:
```bash
git clone --recurse-submodules https://github.com/sanger-pathogens/long-read-ampliseq.git
```

Usage is slightly different:
```bash
nextflow run long-read-ampliseq/main.nf \
--raw_read_dir <directory containing FAST5/POD5 files> \
--reference <reference fasta> \
--primers <fasta containing primers> \
--target_regions_bed <BED file containing target regions> \
--additional_metadata <CSV mapping sample IDs to barcodes> \
--clair3_model <path to Clair3 model> \
-profile standard
```

The standard profile is intended to allow the pipeline to run (with internet access) on the Sanger HPC (farm). It ensures the pipeline can run with the LSF job scheduler and uses singularity images for dependencies management, as well as the latest versions of the pipeline base configuration (from [PaM Info common config file](https://github.com/sanger-pathogens/nextflow-commons/blob/master/configs/nextflow.config)) and Dorado models.

It's best to run the pipeline as a job in the oversubscribed queue i.e. preface the command with this:
```bash
bsub -o output.o -e error.e -q oversubscribed -R "select[mem>4000] rusage[mem=4000]" -M4000
```

Once your job has finished and you're happy with the output, clean up any intermediate files. To do this (assuming no other pipelines are running from the current working directory), run:
```bash
rm -rf work .nextflow*
```

## Support
Please contact PaM Informatics for support through our [helpdesk portal](https://jira.sanger.ac.uk/servicedesk/customer/portal/16) or for external users please reach out by email: pam-informatics@sanger.ac.uk
