Changelog
=========

x.y.z (MM/DD/YYYY)
------------------
- Replaced Tophat2 with HISAT2; removed bbmap.
- Replaced HTSeq with stringtie; consequent renumbering of outputs

  - Added new required config value alignment_option: read_length (and set
    default to 50)

- Replaced Tuxedo/CuffDiff with Ballgown.
- Added a top-level conda environment for watermelon

  - Upgraded Python 3.6.6, Snakemake 5.3.0, pandas (0.23.4)

- Upgraded MultiQC to 1.6 (and adjusted to use conda environment)
- Adjusted config:library_type to accept

  - fr-unstnded
  - unstranded
  - forward_reverse
  - fr-firststrand
  - reverse_forward
  - fr-secondstrand

- Throttled fastqc to avoid Java memory overallocation
- Adjusted watermelon to enable "in-flight" dry-run/dag (executed in the
  directory of a job currently in-progress).


0.3.6 (8/12/2018)
-----------------
- Adjusted how environment modules are versioned and loaded
- Disabled color output to avoid errors on bfx-comp6 (a transient regression, we hope)
- Adjusted watermelon script and rules to work consistently across comp 3,5,6

0.3.5 (8/1/2018)
----------------
- Modularized snakefile by splitting rules into individual files.
- Revised "all" rule to specify the minimal set of outputs.
- Adjusted config to group dirs into single block.
- Fixed bug that caused multiqc to crash when custom alignment dir specified
- Extended version tests to check multiqc installed correctly

0.3.4 (6/12/2018)
-----------------
- Adjusted module files/tests to make compatible with bfx-comp5/6
- Renamed watermelon_rnaseq to watermelon_dependencies
- Adjusted versions of watermelon and watermelon_dependencies modules to match
  Watermelon version number

0.3.3 (12/20/2017)
------------------
- Modified rnaseq.snakefile to wait until all multiqc files are available
  before making the alignment_qc.html

0.3.2 (11/15/2017)
------------------
- Added support for zebrafish (GRCz10)
- Adjusted config validation to fail if test-control comparison values are not distinct
- Added step to create combined gene list summaries in deliverables/
- Modified top 500 gene heatmaps in DESeq2: row scaling, row dendro, and aspect ratio.


0.3.1 (9/25/2017)
-----------------
- Added support for c. elegans (ce10, ce11, WBS235) and GRCh37
- Added support for multiple runs

  - Adjusted watermelon_init to display matrix of sample run files
  - Added validation error where a run or sample has no fastq files
  - Revised how source files are linked during init; hardlinked where
    possible (and symlinked if not)

- Adjusted watermelon to warn and/or gracefully skip DESeq2 if no replicates
  in any phenotype
- Added fastq_screen rule to analyze breakdown of alignments within and
  across species to identify contamination and/or depletion problems
- Adjusted DESeq2 to produce pre and post normalization PCA plots
- Adjusted how genome references are merged with template config to allow for
  nested dicts and also avoid accidentally overwriting default template dicts
- Corrected a bug in DESeq2 MA and volcano plots that incorrectly labeled the
  top 10 diffex genes in PDF output

0.3.0 (7/28/2017)
-----------------
- Added support for paired-end reads
- Transitioned naive alignment QC metrics to MultiQC
- Fixed a bug that caused deseq2 to crash when merging htseq counts files
  with non-standard sample names
- Added memory constraint to watermelon shell script
- Moved references to common location

0.2.5 (6/22/2017)
-----------------
- Fixed DESeq2 bug in correlation plots that crashes watermelon when less
  than 10 samples in constrast
- Simplified DESeq2 plot labels

0.2.4 (5/12/2017)
-----------------
- Adjusted DESeq PCA graphs to include variance percentages in axis labels
- Adjusted DESeq to only consider phenotypes that have replicates (DESeq
  crashes when attempting to normalize phenotypes without replicates)
- Fixed DESeq bug that occurs when phenotype is not compared
- Adjusted tuxedo-cuffdiff gene lists to correctly flip test and controls so
  fold change directions match DESeq2
- Reverted console logging to be verbose thereby avoiding suppression of
  logging under certain error conditions

0.2.3 (5/8/2017)
----------------
- Corrected memory allocation bug in DESeq2/pandoc to prevent DESeq diffex from
  occasionally crashing
- Adjusted watermelon to filter console logging to progress messages
- Speed dry-run mode by skipping follow-on "summary detail" job

0.2.2 (5/3/2017)
----------------
- Corrected bug that crashed pipeline if only one phenotype specified

0.2.1 (4/25/2017)
-----------------
- Corrected the way HTSeq process stranded data
- Adjusted watermelon shell script to always print/log shell commands

0.2 (4/17/2017)
---------------
- Added DESeq2 diffex analysis
  - Adjusted config to include main_factor
  - DESeq2 calling and extensive plots
  - Basic annotation
- Revised and simplified output folders and rule naming
  - tuxedo steps are renumbered
  - config_checksums are hidden
  - log dirs are hidden
  - simplified deliverable rules
- Added diffex comparison gene summaries for tuxedo (cuffdiff) and DESeq2 results
- Improved watermelon launch
  - Config validation check for well-formed (R friendly) phenotype labels,
    and phenotype values
  - Improved handling of locked dir
  - Instead of failing fast on error, watermelon will run the valid remainder of jobs
    (--keep-going)
- Separated legacy dependencies from watermelon dependencies
- Bugfixes:
  - watermelon_init would show confusing result when fastq source was inside working dir
  - Cuffadapt would always run even if cut adapt config params were set to 0
  - HTSeq sometimes failed bc too many threads allocated
  - HTSeq merge sometimes failed due to incomplete files

0.1 (2/10/2017)
---------------
- Initial development release
- Recapitulated legacy functionality in a snakemake implementation
