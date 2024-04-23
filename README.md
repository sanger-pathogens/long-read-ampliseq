# Long read AmpliSeq

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

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
--additional_metadata <CSV of additional metadata>
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
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
