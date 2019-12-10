# Watermelon

An RNA-seq pipeline for the UMich Bioinformatics Core

The official repository is located at [https://github.com/umich-brcf-bioinf/Watermelon](https://github.com/umich-brcf-bioinf/Watermelon)

## Overview

Watermelon interprets fastqs and configuration to produce alignments, QC data, feature counts, and diffex results.
Watermelon uses snakemake, which in turn wraps calls to the various bioinformatic tools
in the pipeline (e.g. cutadapt, fastqc, etc.).

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

The watermelon conda environment has all of the required software to start the pipeline - i.e. to run watermelon_init and snakemake. This environment is unlikely to change frequently, and so will only need to be rebuilt at a maximum of once per release. After the environment is built, it can be activated in any time it needs to be used.
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

It can be deactivated later with

    conda deactivate

## Simple Example - Synthetic reads from this repo

The first step is to run watermelon_init, which will set up the analysis in project directory where you invoke it.

To do this, you will need:

* sample sheet - e.g. config/example_samplesheet.csv
* genome build e.g. GRCh38
* job suffix e.g. _20190821
* paths to sample directories with fastq files

Here's an example (Note: replace /path/to/ with an actual valid path)

    # Create a project directory & navigate there
    mkdir ~/watermelon_example_project
    cd ~/watermelon_example_project
    # Copy the samplesheet from the repo to here
    cp /path/to/Watermelon/config/example_samplesheet.csv .
    # Now run watermelon_init
    watermelon_init --sample_sheet example_samplesheet.csv --genome_build GRCh38 --job_suffix "_$(date +%F)" \
    /path/to/Watermelon/data/sim_reads_human_chr22

Now is a good time to review the output from watermelon_init. It generates the following:

* inputs : Directory with links to the sample dirs of input fastqs
* config.yaml : example configuration file for the pipeline
* watermelon.README : Contains results of input file linking, as well as instructions on modifying the config and running the pipeline

In this example, the config generated should be compatible with the given samples. In a real use-case scenario, the config file should be treated as a template, and edited so that it matches the project plan and dataset that it applies to.

Now you can perform a dry-run, which will validate the config and the workflow.

These same directions should also be found in the watermelon.README file generated by watermelon_init:

    # Start a screen session (for persistence over ssh):
    screen -S watermelon{job_suffix}
    # Activate the conda environment:
    conda activate watermelon
    # To validate the config and check the execution plan:
    snakemake --dryrun --printshellcmds --configfile {config_file} --snakefile Watermelon/rnaseq.snakefile

You should still be in the project directory, and ready to run the pipeline.

To run on bfx-comp5/6 (notice the profile):

    snakemake --use-conda --configfile {config_file} --snakefile Watermelon/rnaseq.snakefile --profile Watermelon/config/profile-comp5-6

To run the pipeline on the GreatLakes compute cluster:

    snakemake --use-conda --configfile {config_file} --snakefile Watermelon/rnaseq.snakefile --profile Watermelon/config/profile-greatlakes

## Further Reading
