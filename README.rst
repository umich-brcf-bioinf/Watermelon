==========
Watermelon
==========

A UM BfxCore RNA-seq pipeline that accepts fastq files and produces diffex results.

The official repository is at:
https://github.com/umich-brcf-bioinf/Watermelon

--------
Overview
--------

Watermelon interprets fastqs and configuration to produce alignments and diffex results.
Watermelon uses Snakemake, which in turn wraps calls to the various bioinformatic tools
in the pipeline (e.g. cutadapt, fastqc, etc.).

There are two main steps:

 * waterlemon-init: Creates template config.yaml file and initializes directory structure.
 * running the pipeline: Using the snakemake utility, config file is validated RNA-seq workflow is evaluated accordingly.


For more information see:

* `QUICKSTART`_ : basic orientation on a simple way to run Watermelon.

* `INSTALL`_ : basic ways to install.

* `DETAILS`_ : details on the output structure, alternative ways to run, and troubleshooting.


--------------------
watermelon_init help
--------------------

::

  $ watermelon_init --help
  usage: watermelon_init [-h] --genome_build
                         [--job_suffix JOB_SUFFIX] --sample_sheet SAMPLE_SHEET
                         source_fastq_dirs [source_fastq_dirs ...]

  Creates template config file and directories for a watermelon rnaseq job.

  positional arguments:
    source_fastq_dirs     One or more paths to run dirs. Each run dir should
                          contain samples dirs which should in turn contain one
                          or more fastq.gz files. The watermelon sample names
                          will be derived from the sample directories.

  optional arguments:
    -h, --help            show this help message and exit
    --genome_build
                          Config template will based on this genome
    --job_suffix JOB_SUFFIX
                          =_08_21 This suffix will be appended to the names of
                          analysis/deliverables dirs and configfile. (Useful if
                          rerunning a job with revised input fastq or revised
                          config; can help differentiate simultaneous watermelon
                          jobs in top/ps.)
    --sample_sheet SAMPLE_SHEET
                          A CSV file containing sample names and
                          phenotype/characteristic information which correspond
                          to these samples. Watermelon_init will verify that
                          sample names in this file match with those found in
                          the input directories.


The output of watermelon init will include an example config file (to be edited according to research project) and a watermelon.README file containing information about how to run the snakemake pipeline.

====

Email bfx-rnaseq-pipeline@umich.edu for support and questions.

UM BRCF Bioinformatics Core

.. _INSTALL: doc/INSTALL.rst
.. _DETAILS: doc/DETAILS.rst
.. _QUICKSTART : doc/QUICKSTART.rst
