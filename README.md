# Long-read Ampliseq

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15181157.svg)](https://doi.org/10.5281/zenodo.15181157)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

A Nextflow pipeline for basecalling, read mapping, QC, variant calling and analysis of nanopore multiplex amplicon data.

![flowchart](images/ampliseq_pipeline.jpg)

For tool versions please see the processes within the 'modules' directory.

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

Should you need to run the pipeline offline, it is best to make use of pre-populated dependency caches. These can be created with any of the supported profiles (e.g. `-profile docker`) by running the pipeline once to completion. You will also need to provide a `--basecall_model_path` (see step 4 [above](#installation))- the laptop profile includes a default local path for this, as well as `--clair3_model` and `--dorado_local_path`.

You can override the default paths using the command line parameters directly when invoking nextflow or by supplying an additional config file in which these parameters are set, using the `-c my_custom.config` nextflow option.

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

Follow steps 5 and 6 [above](#installation) to download a Clair3 model and clone the repo.

Usage is slightly different (you use the standard profile and don't need `--dorado_local_path`):
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
## OS requirements

The pipeline has been tested on multiple operating systems

Linux: Ubuntu 18.04 + 22.04 (Singularity and Docker profiles) Mac OSX (Laptop profile)

## Hardware requirements

This pipeline requires only a standard modern computer with enough RAM/CPU power to support running the standard tools. There are checks within the pipeline to determine max resources avaliable to ensure sensible resource requests.

The most taxing sections of the pipeline are:

- FastQC
- Minimap2
- RAXML-NG
- Clair3

The pipeline will run without GPU support; however, access to GPU hardware (one that is supported by Dorado; see [here}(https://github.com/nanoporetech/dorado?tab=readme-ov-file#platforms)) will result in a much faster runtime.

## Approximate runtime

When executed in **laptop mode**, a full analysis run including basecalling is typically expected to complete in approximately 10–12 hours. However, this estimate is highly sensitive to sequencing depth; larger input datasets can significantly extend basecalling time. To optimize runtime, we recommend adjusting Dorado basecalling accuracy parameters in accordance with your performance requirements (fast, hac, sup)

In **HPC** environments utilizing Singularity containers with GPU acceleration, runtimes are substantially reduced. Under optimal conditions, end-to-end processing has been observed to complete in as little as 40 minutes, and up to 1.5 hours for typical runs. This is however depending on input size and GPU availability.

> **_Note_** On first-time runs, where none of the required tools are cached, the initial setup (including downloading dpendency software and models) introduces an additional overhead of approximately 10 minutes.

## Dependancies

| Tool            | Version | Container                                            |
| --------------- | ------- | ---------------------------------------------------- |
| bcftools        | 1.20    | quay.io/biocontainers/bcftools1.20--h8b25389_0       |
| bedtools        | 2.31.1  | quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1    |
| clair3          | v1.0.9  | hkubal/clair3:v1.0.9                                 |
| pysam           | 0.0.2   | quay.io/sangerpathogens/pysam:0.0.2                  |
| pandas          | 2.2.1   | quay.io/sangerpathogens/pandas:2.2.1                 |
| python_graphics | 1.0.0   | quay.io/sangerpathogens/python_graphics:1.0.0        |
| cutadapt        | 4.7     | quay.io/biocontainers/cutadapt:4.7--py310h4b81fae_1  |
| cuda_dorado     | 0.7.1   | quay.io/sangerpathogens/cuda_dorado:0.7.1 (SEE NOTE) |
| fastqc          | 0.12.1  | quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0      |
| ubuntu          | 20.04   | ubuntu:20.04                                         |
| minimap2        | 2.26    | quay.io/biocontainers/minimap2:2.26--he4a0461_2      |
| multiqc         | 1.22.2  | quay.io/biocontainers/multiqc:1.22.2--pyhdfd78af_0   |
| pod5            | 0.3.6   | quay.io/sangerpathogens/pod5:0.3.6                   |
| pycoqc          | 2.5.2   | quay.io/biocontainers/pycoqc:2.5.2--py_0             |
| samtools        | 1.19.2  | quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1    |
| seqtk           | 1.4     | quay.io/biocontainers/seqtk:1.4--he4a0461_2          |

> **_NOTE_** however we suggest you install your own version of dorado to match your OS

## Support
Please contact PaM Informatics for support through our [helpdesk portal](https://jira.sanger.ac.uk/servicedesk/customer/portal/16) or for external users please reach out by email: pam-informatics@sanger.ac.uk
