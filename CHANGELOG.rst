Changelog
=========

0.x.x (mm/dd/yyyy)
------------------
- Added support for paired-end reads
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
  fold change directios match DESeq2
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
