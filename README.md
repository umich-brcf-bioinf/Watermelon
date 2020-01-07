# Watermelon

[![Build Status](https://travis-ci.com/umich-brcf-bioinf/Watermelon.svg?token=dKagrFefps7qvBArd5yE&branch=master)](https://travis-ci.com/umich-brcf-bioinf/Watermelon)

An RNA-seq pipeline for the UMich Bioinformatics Core

The official repository is located at [https://github.com/umich-brcf-bioinf/Watermelon](https://github.com/umich-brcf-bioinf/Watermelon)

## Overview

Watermelon is an RNAseq pipeline that produces alignments, QC data, feature counts, and differential expression results.
Watermelon uses snakemake, which in turn wraps calls to the various bioinformatic tools in the pipeline.

There are two main steps:

1. Running waterlemon_init: Creates template config.yaml file and initializes input directory structure
2. Running the pipeline: Using the snakemake utility, config file is validated and RNA-seq workflow is evaluated accordingly.

## Getting Started - Installation

UMich Bioinformatics Core and Advanced Genomics Core can use the installation located at:

    /nfs/med-bfx-common/pipelines/Watermelon/Watermelon-seedless/

This is kept up to date with the master branch of this repo, so installation is already complete

Alternatively, you can clone this repository

    # To install in your home directory
    cd ~/
    git clone https://github.com/umich-brcf-bioinf/Watermelon

In order to have access to watermelon_init (which will set up an example config and prepare for the pipeline to be run),
you must make sure that the watermelon/bin location is in your PATH environment variable.
To achieve this one-time, execute the following. For permanent effect, place the same command in your ~/.bash_profile or ~/.bashrc

    # If you're UMich BfxCore or Advanced Genomics Core:
    PATH=$PATH:/nfs/med-bfx-common/pipelines/Watermelon/Watermelon-seedless/bin/
    # Or if you've cloned to your home directory
    PATH=$PATH:/home/$USER/Watermelon/bin/

## Getting Started - Conda environment

The watermelon conda environment has all of the required software to start the pipeline. This environment is unlikely to change frequently, and so will only need to be rebuilt at a maximum of once per release. After the environment is built, it can be activated whenever you want to run the pipeline, and later deactivated when it is no longer needed.
Note: If you don't already have anaconda3/miniconda3, then I'd recommend installing miniconda3 - [Link to conda installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

There is a .yaml file in this repo which can be used to build the watermelon conda environment

    # Navigate to Watermelon (BfxCore/MAGC)
    cd /nfs/med-bfx-common/pipelines/Watermelon/Watermelon-seedless/
    # Alternatively, if you cloned it to your home dir, navigate there
    cd ~/Watermelon
    # Create the conda environment from the .yaml file
    conda env create -f envs/watermelon.yaml

After the conda environment is created, you can activate it with

    conda activate watermelon

This environment provides the basic software requirements to run watermelon_init and snakemake. It can be deactivated later with

    conda deactivate

## Walkthrough - Example w/ simulated reads

This repository contains a set of simulated paired-end read data which will be used for this example. The provided example configuration file and samplesheet used here are set up to run this on this example without any modification, for illustrative purposes. For more examples, see the "Further Reading" section below.

The first step is to run watermelon_init, which will set up the analysis in the project directory where you invoke it.

To do this, you will need:

* sample sheet
  * e.g. config/example_samplesheet.csv
  * must have column labeled `sample` (all lowercase) with sample IDs as the values
  * optional additional columns with phenotype or other sample info for differential expression analysis
* genome build
  * e.g. GRCh38
  * Available genomes `GRCh38, GRCh37, GRCm38, NCBIM37, Rnor_6.0, GRCz11, BDGP6, WBcel235`
  * These are ENSEMBL reference IDs. [Table of equivalent UCSC  & ENSEMBL IDs](doc/Equivalence_UCSC_ENSEMBL.md)
  * There are also the special genome build options `Other` & `TestData`.
    * The link to alternative reference example below outlines using `Other` for customized references.
    * The current example uses the `TestData` genome build, which is a human chr22 reference included with this repository.
* job suffix
  * e.g. `_20190821`
  * Suffix will be applied to config name and analysis results dir
* paths to sample directories with fastq files
  * e.g. `/path/to/Watermelon/data/sim_reads_human_chr22`
  * Path given on command line must point to a directory containing subdirectories with names matching samples in samplesheet
  * Fastq files for a sample will reside in the appropriately named subdirectory
  * Fastq files must have `_R1.fastq[.gz]` or `_R2.fastq[.gz]` in their filenames


Here's an example of running waterlemon_init (Note: replace /path/to/Watermelon in the following code block with valid paths)

    # Create a project directory & navigate there
    mkdir ~/watermelon_example_project
    cd ~/watermelon_example_project
    # Copy the samplesheet from the repo to here
    cp /path/to/Watermelon/config/example_samplesheet.csv .
    # Now run watermelon_init
    watermelon_init --sample_sheet example_samplesheet.csv --genome_build GRCh38 --job_suffix _20190821 \
    /path/to/Watermelon/data/sim_reads_human_chr22

Now is a good time to review the output from watermelon_init. It generates the following:

* inputs : Directory with links to the sample dirs of input fastqs
* config.yaml : example configuration file for the pipeline
* watermelon.README : Contains results of input file linking, as well as instructions on modifying the config and running the pipeline

In this example, the config generated should be compatible with the given samples. In a real use-case scenario, the config file should be treated as a template, and edited so that it matches the project plan and dataset that it applies to.

Now you can perform a dry-run, which will validate the config and the workflow.

These same directions should also be found in the watermelon.README file generated by watermelon_init:

    # Start a screen session (for persistence over ssh):
    screen -S watermelon_20190821
    # Activate the conda environment:
    conda activate watermelon
    # Dry-run to validate the config and check the execution plan:
    snakemake --dryrun --printshellcmds --configfile config_20190821.yaml --snakefile Watermelon/rnaseq.snakefile

You should still be in the project directory, and ready to run the pipeline.

To run on bfx-comp5/6 (notice the profile):

    snakemake --use-conda --configfile config_20190821.yaml --snakefile Watermelon/rnaseq.snakefile --profile Watermelon/config/profile-comp5-6

Similarly, to run the pipeline on the GreatLakes compute cluster:

    snakemake --use-conda --configfile config_20190821.yaml --snakefile Watermelon/rnaseq.snakefile --profile Watermelon/config/profile-greatlakes

## Further Reading

* [Troubleshooting](doc/troubleshooting.md)
* [Pipeline rulegraph](doc/rulegraph.svg)
* [Example - MAGC data](doc/example_magc_data.md)
* [Example - Alignment, feature count, & QC Only](doc/example_align_qc_only.md)
* [Example - Alternative references](doc/example_alt_refs.md)
* [Example - Generating an annotation TSV file](doc/generating_annotation_tsv.md)
* [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/)
