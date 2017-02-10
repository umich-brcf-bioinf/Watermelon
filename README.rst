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
Watermelon wraps Snakemake, which in turn wraps calls to the various bioinformatic tools
in the pipeline (e.g. cutadapt, tophat, etc.).

There are two executables:

# waterlemon-init: Creates template config.yaml file and initializes directory structure.
# watermelon: Interprets the config file and executes the RNA-seq workflow.
              Facilitates validation, logging, and sets basic default options.

Watermelon workflow consists of two main steps:

   * sequencing [fastqs] -> watermelon-init -> config.yaml
   * [manual edits to config file] -> watermelon -> results


For more information see:

* `QUICKSTART`_ : basic orientation on a simple way to run Watermelon.

* `INSTALL`_ : basic ways to install.

* `DETAILS`_ : details on the output structure, alternative ways to run, and troubleshooting.


--------------------
watermelon-init help
--------------------

::

  $ watermelon-init --help
   usage: watermelon-init [-h] --genome_build {hg19,mm10,rn5}
                          [--job_suffix JOB_SUFFIX]
                          source_fastq_dir
   
   Creates template config file and directories for a watermelon rnaseq job.
   
   positional arguments:
     source_fastq_dir      Path to dir which contains samples dirs (each sample
                           dir contains one or more fastq.gz files). The sample
                           dir names are used to initialize the config template.
                           If the source_fastq_dir is outside the working dir,
                           init will make local inputs dir containing symlinks to
                           the original files.
   
   optional arguments:
     -h, --help            show this help message and exit
     --genome_build {hg19,mm10,rn5}
                           Config template will based on this genome
     --job_suffix JOB_SUFFIX
                           =_02_10 This suffix will be appended to the names of
                           analysis/deliverables dirs and configfile. (Useful if
                           rerunning a job with revised input fastq or revised
                           config; can help differentiate simultaneous watermelon
                           jobs in top/ps.)
   

--------------------
watermelon help
--------------------

::

  $ watermelon-init --help
   Usage: watermelon [options] -c {config_file}
   Example: watermelon -c ~/my_config.yaml
   
   Executes the RNA-seq workflow by wrapping the Snakemake call, specifically setting common,
   useful defaults and passing through any unrecognized command line options through to
   snakemake. Captures all outputs in a log and times overall execution.
   
   Options:
       -c, --configfile : [config.yaml] snakemake config file
       --cores [N]      : =40, use at most N cores in parallel
       --dag            : create DAG files (in working dir) to visualize execution plan
       -n, --dryrun     : show plan without executing anything
       --help           : shows this message


====

Email bfx-rnaseq-pipeline@umich.edu for support and questions.

UM BRCF Bioinformatics Core

.. _INSTALL: doc/INSTALL.rst
.. _DETAILS: doc/DETAILS.rst
.. _QUICKSTART : doc/QUICKSTART.rst