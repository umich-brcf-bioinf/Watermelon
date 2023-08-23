# Watermelon

An RNA-seq pipeline for the UMich Bioinformatics Core

The official repository is located at [https://github.com/umich-brcf-bioinf/Watermelon](https://github.com/umich-brcf-bioinf/Watermelon)

* [Overview](#overview)
* [Getting Started - Installation](#getting-started---installation)
* [Getting Started - Conda Environment](#getting-started---conda-environment)
* [Walkthrough - Alignment and QC Example](#walkthrough---alignment-and-qc-example)
* [Walkthrough - Differential Expression Example](#walkthrough---differential-expression-example)
* [Further Reading](#further-reading)

## Overview

Watermelon is an RNAseq pipeline that produces alignments, QC data, feature counts, and differential expression results.
Watermelon uses snakemake, which in turn wraps calls to the various bioinformatic tools in the pipeline.

There are two main steps:

1. Running waterlemon_init: Creates template config.yaml file
2. Running the pipeline: Using the snakemake utility, an RNA-seq workflow (either alignment + QC workflow or differential expression workflow) is evaluated accordingly.

## Getting Started - Installation

UMich Bioinformatics Core and Advanced Genomics Core can use the installation located at:

    /nfs/turbo/umms-brcfpipeline/pipelines/Watermelon/

This is kept up to date with the master branch of this repo, so installation is already complete

Alternatively, you can clone this repository

    # To install in your home directory
    cd ~/
    git clone https://github.com/umich-brcf-bioinf/Watermelon

In order to have access to watermelon_init (which will set up an example config and prepare for the pipeline to be run),
you must make sure that the watermelon location is in your PATH environment variable.
To achieve this one-time, execute the following. For permanent effect, place the same command in your ~/.bash_profile or ~/.bashrc

    # If you're UMich BfxCore or Advanced Genomics Core:
    PATH=$PATH:/nfs/turbo/umms-brcfpipeline/pipelines/Watermelon/
    # Or if you've cloned to your home directory
    PATH=$PATH:/home/$USER/Watermelon/

## Getting Started - Conda Environment

The watermelon conda environment has all of the required software to start the pipeline. This environment is unlikely to change frequently, and so will only need to be rebuilt at a maximum of once per release. After the environment is built, it can be activated whenever you want to run the pipeline, and later deactivated when it is no longer needed.
Note: If you don't already have anaconda3/miniconda3, then I'd recommend installing miniconda3 - [Link to conda installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

There is a .yaml file in this repo which can be used to build the watermelon conda environment

    # Navigate to Watermelon (BfxCore/MAGC)
    cd /nfs/turbo/umms-brcfpipeline/pipelines/Watermelon/
    # Alternatively, if you cloned it to your home dir, navigate there
    cd ~/Watermelon
    # Create the conda environment from the .yaml file
    conda env create -f envs/watermelon.yaml

After the conda environment is created, you can activate it with

    conda activate watermelon

This environment provides the basic software requirements to run watermelon_init and snakemake. It can be deactivated later with

    conda deactivate

## Walkthrough - Alignment and QC Example

This repository contains a set of simulated paired-end read data which will be used for this example. The provided example configuration file and samplesheet used here are set up to run this on this example without any modification, for illustrative purposes. For more examples, see the "Further Reading" section below.

The first step is to run watermelon_init, which will set up the analysis in the project directory where you invoke it.

For this example, you will need to provide:

* genome build
  * e.g. GRCh38
  * Available genomes `GRCh38, GRCh37, GRCm38, NCBIM37, Rnor_6.0, GRCz11, BDGP6, WBcel235`
  * These are ENSEMBL reference IDs. [Table of equivalent UCSC  & ENSEMBL IDs](doc/Equivalence_UCSC_ENSEMBL.md)
  * There are also the special genome build options `Other` & `TestData`
    * The link to alternative reference example below outlines using `Other` for customized references
    * The current example uses the `TestData` genome build, which is a human chr22 reference included with this repository.
* project ID
  * e.g. given `20190821`, a `_20190821` suffix will be applied to config name and analysis results dir
* config type
  * Either `align_qc` or `diffex`. Will produce a different config file depending on the given type
* path to at least one sample directory with fastq files
  * e.g. `/path/to/Watermelon/data/sim_reads_human_chr22`
  * Path given on command line must point to a directory containing fastqs. Watermelon_init will attempt to create a shell glob pattern for each of the samples and place this into the generated samplesheet
  * Fastq files must have `_R1` or `_R2` in their filenames to distinguish between read 1 and read 2. Good filename endings would be `_R1.fastq.gz` and `_R2.fastq.gz`


Here's an example of running waterlemon_init (Note: replace /path/to/Watermelon in the following code block with valid paths)

    # Start a screen session (for persistence over ssh):
    screen -S watermelon_20190821
    # Activate the conda environment:
    conda activate watermelon
    # Create a project directory & navigate there
    mkdir ~/watermelon_example_project
    cd ~/watermelon_example_project
    # Now run watermelon_init
    watermelon_init.py --genome_build TestData --project_id 20190821 --type align_qc --input_run_dirs /path/to/Watermelon/data/sim_reads_human_chr22

Now is a good time to review the output from watermelon_init. It generates the following:

* config_20190821.yaml : example configuration file for the pipeline
* samplesheet.csv : A CSV samplesheet (auto-generated from the input_run_dirs contents) containing a `sample` column and a `run_dir` column to be used by the pipeline
* Watermelon : Directory containing a copy of the pipeline code


Now you can perform a dry-run, which will validate that the config is compatible with the workflow.

    # Singularity must be available to snakemake, for environment management under the hood
    module load singularity
    # Dry-run to validate the config and check the execution plan:
    snakemake --dryrun --printshellcmds --configfile config_20190821.yaml --snakefile Watermelon/align_qc.smk

You should still be in the project directory, and ready to run the pipeline.

To run on bfx-comp5/6 (notice the profile):

    snakemake --configfile config_20190821.yaml --snakefile Watermelon/align_qc.smk --profile Watermelon/config/profile-comp5-6

Similarly, to run the pipeline on the GreatLakes compute cluster:

    # Singularity must be available to snakemake, for environment management under the hood
    module load singularity
    snakemake --configfile config_20190821.yaml --snakefile Watermelon/align_qc.smk --profile Watermelon/config/profile-greatlakes

## Walkthrough - Differential Expression Example

A differential expression workflow starts out the same, but differs in execution:

1. Use watermelon_init to create a config file for the differential expression analysis
2. Inspect the generated config file, modify differential expression details based on experimental design
3. Run the deseq2_analysis.R script using the WAT_diffex singularity image, which provides a standardized compute environment

When running watermelon_init, as before you'll provide a genome build, a project ID, and a config type. Since it will be the `diffex` type, we will also give it a count matrix and a samplesheet. I will take the same samplesheet that was generated earlier, and add an additional column `treatment` with half of the samples labeled `control` and the other half `drug`. In this example, the config generated will be compatible with the given samples. In a real use-case scenario, the config file should be treated as a template, and edited so that it matches the project plan and dataset that it applies to.

Before running watermelon_init, I'll grab the count matrix (excluding annotation columns).

    cut -f1,5- analysis_20190821/deliverables/counts/gene_expected_count.annot.txt > counts_ready.txt

Count matrix details:
- Text file of tab-separated values
- Row for each gene, matching the IDs of the GTF used in genome_build
- Column for each sample with count values, with column names matching the samples in the samplesheet

Note:
It is important that the count matrix only contains the above information. If annotation columns are present (such as description), these must be removed before running the pipeline

Running watermelon_init for `diffex` type:

    watermelon_init.py --genome_build TestData --project_id 20190821d --type diffex --sample_sheet samplesheet.csv --count_matrix counts_ready.txt

Notes: It will prompt about overwriting Watermelon. In this example, it's inconsequential. In real use-cases, use your judgement and move/rename things if overwrite isn't desired. Also, you might've noticed the slightly altered project_id in this example (suffixed with d) to prevent overwriting the existing config and/or mixing the results from the workflows.

After running this, you'll have:

* config_20190821d.yaml : a `diffex` type configuration file for the differential expression workflow. In this example, it needs no modification to work with the differential expression workflow. In a real use-case scenario, you will modify this config file according to the project plan.
* Watermelon : Directory containing a copy of the pipeline code

For this example we should modify our samplesheet that we auto-generated earlier, so that it matches our example analysis. In a real example, you'd likely have your samplesheet already set up with phenotypes and conditions, and you would then modify the config to match. Here, our example analysis compares two `treatment`s, `control` and `drug`. So let's add a column to our samplesheet with this information. We could remove or exclude the input_glob column, since it's not used in the diffex analysis, but we'll leave it in place now because that's slightly easier.

    sample,input_glob,treatment
    sample_01,...,control
    sample_02,...,control
    sample_03,...,control
    sample_04,...,control
    sample_05,...,drug
    sample_06,...,drug
    sample_07,...,drug
    sample_08,...,drug


Now we can run the deseq2_analysis.R script. Since we're using a singularity image, we'll first set an environment variable that will define a read-only bind mount for the annotation data.

    # Define the read-only bind mount to our reference data
    export SINGULARITY_BIND="/nfs/turbo/umms-brcfpipeline/references:/nfs/turbo/umms-brcfpipeline/references:ro"
    # Run the deseq2_analysis.R script, providing environment with the singularity image
    singularity exec docker://umichbfxcore/wat_diffex:0.6.0 Rscript Watermelon/scripts/deseq2_analysis.R --configfile config_20190821d.yaml --markdownfile Watermelon/report/report_diffex.Rmd 2>&1 | tee deseq2_$(date +%FT%H%M%S).log

When this runs, by default it will run the DESeq2 analysis and save the outputs, including a `.Rdata` file, and then it will knit a report using these outputs. 


> Note: If it is desired to run just the analysis portion of the script or just the reporting portion of the script, `deseq2_analysis.R` has command line flags to allow this behavior, `--no_knit` and `--no_analysis`, respectively.

> Note: It is also possible to use the singularity image with an RStudio session, similar to how it's described on our [single cell analysis environment page](https://github.com/umich-brcf-bioinf/single_cell_env), if an interactive session is desired.

Next we will edit (if desired) and finalize the report. Inspect `report_draft.html`, and make any desired edits to the `report_draft.md` that you need. After you're satisfied with the edits, knit this `report_draft.md` into our final report.

    # Run the deseq2_analysis.R script again, this time with
    # our report_draft.md as the --markdownfile
    # and with the flag --report_finalize
    #
    singularity exec docker://umichbfxcore/wat_diffex:0.6.0 Rscript Watermelon/scripts/deseq2_analysis.R --configfile config_20190821d.yaml --markdownfile analysis_20190821d/report/report_draft.md --report_finalize


Finally, we will move the deliverables into their desired location. Within the DESeq2 analysis script, we've tracked all deliverables and created a file - `deliverables_list.txt` - that can be used directly with `rsync` to transfer results to a desired location.

    rsync --files-from deliverables_list.txt . analysis_20190821d/deliverables

Files can be delivered anywhere - just replace final argument with your destination (no trailing slash!)

    # rsync --files-from deliverables_list.txt . <dest>

> Note: The `deliverables_list.txt` file has `/./` that are inserted into the filepaths by the analysis script when it is written. This is a syntax recognized and used by rsync, where directories are re-created after the `/./` but not before, allowing the desired output structure to be placed anywhere.



## Further Reading

* [Alignment & QC Report generation](doc/report_generation.md)
* [Troubleshooting](doc/troubleshooting.md)
* [Example - Multiple sequencing runs](doc/example_multiple_runs.md)
* [Example - Alternative references](doc/example_alt_refs.md)
* [Example - Generating an annotation TSV file](doc/generating_annotation_tsv.md)
* [Notes - Generating or updating the references](doc/creating_updated_refs.md)
* [Notes - Git hooks for automated tests](doc/git_hooks_auto_tests.md)
* [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/)
