==========
Watermelon
==========

An bioinformatics pipeline to analyze RNASeq data.

--------
Overview
--------

Accepts a set of fastq.gz files and trims, QC, align, calls diffex genes, and generates reports.


----------
Quickstart
----------

TODO

----------
Files/dirs
----------
* config/ 
  sample config.yamls for common experiment types
* custom_lib
  3rd party scripts that required customization to run correctly/deterministically
* doc
  TODO
* modulefiles
  environment module definitions relevant to the pipeline
* rnaseq.snakefile
  top-level snakefile for rnaseq processing
* scripts
  custom scripts referenced by snakefile above
