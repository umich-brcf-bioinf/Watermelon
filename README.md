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
  * There are also the special genome build options `Other` & `TestData`.
    * The link to alternative reference example below outlines using `Other` for customized references.
    * The current example uses the `TestData` genome build, which is a human chr22 reference included with this repository.
* project ID
  * e.g. given `20190821`, a `_20190821` suffix will be applied to config name and analysis results dir
* config type
  * Either `align_qc` or `diffex`. Will produce a different config file depending on the given type
* paths to sample directories with fastq files
  * e.g. `/path/to/Watermelon/data/sim_reads_human_chr22`
  * Path given on command line must point to a directory containing subdirectories with names matching samples in samplesheet
  * Fastq files for a sample will reside in the appropriately named subdirectory
  * Fastq files must have `_R1.fastq[.gz]` or `_R2.fastq[.gz]` in their filenames


Here's an example of running waterlemon_init (Note: replace /path/to/Watermelon in the following code block with valid paths)

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

    # Start a screen session (for persistence over ssh):
    screen -S watermelon_20190821
    # Activate the conda environment:
    conda activate watermelon
    # Dry-run to validate the config and check the execution plan:
    snakemake --dryrun --printshellcmds --configfile config_20190821.yaml --snakefile Watermelon/align_qc.smk

You should still be in the project directory, and ready to run the pipeline.

To run on bfx-comp5/6 (notice the profile):

    # Singularity must be available to snakemake, for environment management under the hood
    module load singularity
    snakemake --configfile config_20190821.yaml --snakefile Watermelon/align_qc.smk --profile Watermelon/config/profile-comp5-6

Similarly, to run the pipeline on the GreatLakes compute cluster:

    # Singularity must be available to snakemake, for environment management under the hood
    module load singularity
    snakemake --configfile config_20190821.yaml --snakefile Watermelon/align_qc.smk --profile Watermelon/config/profile-greatlakes

## Walkthrough - Differential Expression Example

With a very similar process, we can run the differential expression workflow:

1. Use watermelon_init to create a file for the differential expression workflow
2. Inspect the generated config file, modify differential expression details based on experimental design
3. Run the differential expression snakemake workflow with the newly created config file

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

* config_20190821d.yaml : a `diffex` type configuration file for the differential expression workflow. In this example, it needs no modification to work with the differential expression workflow.
* Watermelon : Directory containing a copy of the pipeline code


Now we can run a dry-run similar to above:

    # Assuming still in screen session with watermelon conda env activated
    snakemake --dryrun --printshellcmds --configfile config_20190821d.yaml --snakefile Watermelon/deseq2.smk

If that works fine, then the pipeline can be run. Example of running this on Great Lakes cluster:

    module load singularity
    snakemake --configfile config_20190821d.yaml --snakefile Watermelon/deseq2.smk --profile Watermelon/config/profile-greatlakes



## Further Reading

* [Report generation](doc/report_generation.md)
* [Troubleshooting](doc/troubleshooting.md)
* [Example - Multiple sequencing runs](doc/example_magc_data.md)
* [Example - Alternative references](doc/example_alt_refs.md)
* [Example - Generating an annotation TSV file](doc/generating_annotation_tsv.md)
* [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/)
