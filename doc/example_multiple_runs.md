# Example - Running with multiple sequencing runs from UMich Advanced Genomics Core

The walkthrough in the README provided an example that should work out-of-the-box, where things were simplified for illustrative purposes. The current example goes a bit further, highlighting a few more details which might be useful when running with multiple sequencing runs.

## Before starting:

Just as in the readme walkthrough, you'll run watermelon_init to set up the project directory. Also like before, the watermelon conda environment must be activated in order to run watermelon_init, and singularity must be available.

## Multiple sequencing runs, different biological samples

This first scenario is the simplest - when you want to combine sequencing runs in a workflow, but the samples included in each run are distinct. In this case, you can provide watermelon_init with multiple input_run_dirs, and each sample will be mapped to its appropriate input directory in the samplesheet.

    watermelon_init.py --genome_build GRCh38 --project_id 20190822 --type align_qc --input_run_dirs /path/to/NovaA_1234 /path/to/NovaA_1256

With this command, watermelon_init will create a samplesheet that looks something like the following:

    sample,input_glob
    1111-AB-1,/path/to/NovaA_1234/fastqs_1111-AB/1111-AB-1_*.fastq.gz
    1111-AB-2,/path/to/NovaA_1234/fastqs_1111-AB/1111-AB-2_*.fastq.gz
    1111-AB-3,/path/to/NovaA_1234/fastqs_1111-AB/1111-AB-3_*.fastq.gz
    1111-AB-4,/path/to/NovaA_1234/fastqs_1111-AB/1111-AB-4_*.fastq.gz
    1111-AB-5,/path/to/NovaA_1234/fastqs_1111-AB/1111-AB-5_*.fastq.gz
    1111-AB-6,/path/to/NovaA_1234/fastqs_1111-AB/1111-AB-6_*.fastq.gz
    1111-AB-7,/path/to/NovaA_1234/fastqs_1111-AB/1111-AB-7_*.fastq.gz
    1111-AB-8,/path/to/NovaA_1234/fastqs_1111-AB/1111-AB-8_*.fastq.gz
    1111-AB-21,/path/to/NovaA_1256/fastqs_1111-AB/1111-AB-21_*.fastq.gz
    1111-AB-22,/path/to/NovaA_1256/fastqs_1111-AB/1111-AB-22_*.fastq.gz
    1111-AB-23,/path/to/NovaA_1256/fastqs_1111-AB/1111-AB-23_*.fastq.gz
    1111-AB-24,/path/to/NovaA_1256/fastqs_1111-AB/1111-AB-24_*.fastq.gz
    1111-AB-25,/path/to/NovaA_1256/fastqs_1111-AB/1111-AB-25_*.fastq.gz
    1111-AB-26,/path/to/NovaA_1256/fastqs_1111-AB/1111-AB-26_*.fastq.gz
    1111-AB-27,/path/to/NovaA_1256/fastqs_1111-AB/1111-AB-27_*.fastq.gz
    1111-AB-28,/path/to/NovaA_1256/fastqs_1111-AB/1111-AB-28_*.fastq.gz

In this process, watermelon_init has made some minimal inferences to determine an appropriate shell glob for each sample. In the above example, the samplesheet that was created by watermelon_init should have the right contents so that it can be directly used by the pipeline. When the pipeline runs, these shell globs will expand into the list of fastq files that are used as inputs for a given sample. As noted in the README, the fastq files must contain a `_R1` or `_R2` in their filenames, because while running, the pipeline uses these identifiers to treat the files appropriately. For example, R1s will only be concatenated with R1s, likewise R2s only with R2s, several steps will run in 'paired-end mode' if these are detected, and so forth.

## Multiple sequencing runs, same biological samples

The next scenario is a just a bit more involved in setting up. Here we'll have a set of identical biological samples that were sequenced across two separate runs. This is sometimes done when more sequencing depth is desired. In this case, we want to merge the samples across the runs, so that we can effectively achieve greater sequencing depth by concatenating the samples' fastq files before alignment. I'll demonstrate how to handle this.

In this scenario, NovaA_3456 contains paired-end sequencing data for several samples, and NovaA_5567 contains a subset of those same samples that were resubmitted for sequencing at a later date.

    contents of NovaA_3456/fastqs_3333-CD/
    ---------------------
    NovaA_3456/fastqs_3333-CD/
    ├── 3333-CD-1_CTACGACA-TTGGACTC_S1_R1.fastq.gz
    ├── 3333-CD-1_CTACGACA-TTGGACTC_S1_R2.fastq.gz
    ├── 3333-CD-2_GGACTTGG-CGTCTGCG_S2_R1.fastq.gz
    ├── 3333-CD-2_GGACTTGG-CGTCTGCG_S2_R2.fastq.gz
    ├── 3333-CD-3_TAAGTGGT-GGCTTAAG_S3_R1.fastq.gz
    ├── 3333-CD-3_TAAGTGGT-GGCTTAAG_S3_R2.fastq.gz
    ├── 3333-CD-4_ATATGGAT-TAATACAG_S4_R1.fastq.gz
    ├── 3333-CD-4_ATATGGAT-TAATACAG_S4_R2.fastq.gz
    ├── 3333-CD-5_TCTCTACT-GAACCGCG_S5_R1.fastq.gz
    ├── 3333-CD-5_TCTCTACT-GAACCGCG_S5_R2.fastq.gz
    ├── 3333-CD-6_CCAAGTCT-TCATCCTT_S6_R1.fastq.gz
    ├── 3333-CD-6_CCAAGTCT-TCATCCTT_S6_R2.fastq.gz
    ├── 3333-CD-7_CGGCGTGA-ACAGGCGC_S7_R1.fastq.gz
    ├── 3333-CD-7_CGGCGTGA-ACAGGCGC_S7_R2.fastq.gz
    ├── 3333-CD-8_ATGTAAGT-CATAGAGT_S8_R1.fastq.gz
    └── 3333-CD-8_ATGTAAGT-CATAGAGT_S8_R2.fastq.gz

    contents of NovaA_5567/fastqs_3333-CD/
    ---------------------
    NovaA_5567/fastqs_3333-CD/
    ├── 3333-CD-1_CTACGACA-TTGGACTC_S1_R1.fastq.gz
    ├── 3333-CD-1_CTACGACA-TTGGACTC_S1_R2.fastq.gz
    ├── 3333-CD-2_GGACTTGG-CGTCTGCG_S2_R1.fastq.gz
    ├── 3333-CD-2_GGACTTGG-CGTCTGCG_S2_R2.fastq.gz
    ├── 3333-CD-3_TAAGTGGT-GGCTTAAG_S3_R1.fastq.gz
    ├── 3333-CD-3_TAAGTGGT-GGCTTAAG_S3_R2.fastq.gz
    ├── 3333-CD-4_ATATGGAT-TAATACAG_S4_R1.fastq.gz
    └── 3333-CD-4_ATATGGAT-TAATACAG_S4_R2.fastq.gz

In a general sense, we can begin by running watermelon_init with multiple `input_run_dirs` as we did in the above scenario. However, there are a couple of things to be mindful of before and after watermelon_init runs, which are the following:

* It is imperative that these input directories share a common base path.
    * A recommendation is to place both in the project directory in parallel
    * Valid: `/one/path/to/NovaA_3456` & `/one/path/to/NovaA_5567`
    * Invalid: `/this/path/for/NovaA_3456` & `/another/different/path/to/NovaA_5567`
* The filenames for the same biological samples must have the same prefix. For example, sample 1 in both runs should start with `3333-CD-1_`
* After running watermelon_init, we should modify the samplesheet shell globs so that they expand to the proper set of fastq files for each sample

Knowing this, let's run watermelon_init and then continue from there:

    #In the project directory, run watermelon_init
    watermelon_init.py --genome_build GRCh38 --project_id 20190822 --type align_qc --input_run_dirs NovaA_3456/ NovaA_5567/

This time, watermelon_init will produce a samplesheet that we must modify before we continue:

    sample,input_glob
    3333-CD-1,/path/to/project/NovaA_3456/fastqs_3333-CD/3333-CD-1_*.fastq.gz
    3333-CD-2,/path/to/project/NovaA_3456/fastqs_3333-CD/3333-CD-2_*.fastq.gz
    3333-CD-3,/path/to/project/NovaA_3456/fastqs_3333-CD/3333-CD-3_*.fastq.gz
    3333-CD-4,/path/to/project/NovaA_3456/fastqs_3333-CD/3333-CD-4_*.fastq.gz
    3333-CD-5,/path/to/project/NovaA_3456/fastqs_3333-CD/3333-CD-5_*.fastq.gz
    3333-CD-6,/path/to/project/NovaA_3456/fastqs_3333-CD/3333-CD-6_*.fastq.gz
    3333-CD-7,/path/to/project/NovaA_3456/fastqs_3333-CD/3333-CD-7_*.fastq.gz
    3333-CD-8,/path/to/project/NovaA_3456/fastqs_3333-CD/3333-CD-8_*.fastq.gz

You should notice that these shell globs will not produce our desired set of fastq files when expanded; they only include the fastqs in the first run directory - `NovaA_3456`. However, since we followed the recommendations above, we can modify the shell globs in the samplesheet so that they expand to all of the fastqs that we want:

    3333-CD-1,/path/to/project/NovaA_*/fastqs_3333-CD/3333-CD-1_*.fastq.gz
    3333-CD-2,/path/to/project/NovaA_*/fastqs_3333-CD/3333-CD-2_*.fastq.gz
    3333-CD-3,/path/to/project/NovaA_*/fastqs_3333-CD/3333-CD-3_*.fastq.gz
    3333-CD-4,/path/to/project/NovaA_*/fastqs_3333-CD/3333-CD-4_*.fastq.gz
    3333-CD-5,/path/to/project/NovaA_*/fastqs_3333-CD/3333-CD-5_*.fastq.gz
    3333-CD-6,/path/to/project/NovaA_*/fastqs_3333-CD/3333-CD-6_*.fastq.gz
    3333-CD-7,/path/to/project/NovaA_*/fastqs_3333-CD/3333-CD-7_*.fastq.gz
    3333-CD-8,/path/to/project/NovaA_*/fastqs_3333-CD/3333-CD-8_*.fastq.gz

After making this small change to the samplesheet and verifying that the shell globs expand in the way that we desire, we can continue running the pipeline in the same way we usually do, starting with a dry-run, and then running the pipeline:

    # Dry-run to validate the config and check the execution plan:
    snakemake --configfile config_20190822.yaml --snakefile Watermelon/align_qc.smk -n -p
    # Run the pipeline
    snakemake --configfile config_20190822.yaml --snakefile Watermelon/align_qc.smk --profile Watermelon/config/profile-comp5-6
