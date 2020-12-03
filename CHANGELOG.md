# Changelog

1.6.0 (12/03/2020)
------------------
- Added standalone methods doc to report section & deliverables
- Bugfix - MultiQC silent out-of-memory on cluster
- Bugfix - Hard-linked deliverable intermediates, trimmed reads & aligned bams, are problematic for snakemake

1.5.3 (11/11/2020)
------------------
- Added trimmed fastqs and sorted bams to deliverables
- Updated job resource configurations
- Config validation now accepts 'y' or 'n' answers

1.5.2 (10/23/2020)
------------------
- Reporting improvements from AGC suggestions
    - Wording changes
    - Flexibility (authorship, configurable wording, etc)
    - Some automatic options w/ --AGC flag of init
- Resouce adjustments on several rules
- Consolidated log locations - in `job_logs`
- Documentation for generating reports (for analysts)
- Table of contents in Readme
- Bugfix - samplesheets w/ numeric entries caused uninformative errors
- Bugfix - fastq_screen deliverables rule

1.5.1 (09/08/2020)
------------------
- watermelon_init fix: only use turbo space for default refs now
- Improvements for running on GreatLakes
    - Increased time for trimming
    - Fix for early cutadapt termination recovery

1.5.0 (08/13/2020)
------------------
- Automatic report generation
	- Different report generated for each scenario (align/QC + diffex, align/QC only, diffex only)
	- Numerous tables & figures included in report
	- Links to colocated multiQC report
- Sliced Watermelon goals actualized
	- To run alignment + QC + diffex, no change (default)
	- To run alignment + QC only, either delete entire `diffex` section from config OR add `--config diffex=''` to snakemake command.
	- To run diffex starting from a count matrix, either `count_matrix: /path/to/count_matrix.tsv` must be in config or add `--config count_matrix=/path/to/count_matrix.tsv` to snakemake command.
    - Note: Supplying watermelon_init with `--count_matrix /path/to/count_matrix.tsv` instead of a list of fastq input directories creates an appropriate config containing this key-value pair.
- Genome references moved to turbo location (will be used in coordination w/ AGC)
	- Adjustments to singularity arguments to accommodate this
- Improvements to config validator
	- Use jsonschema module instead of snakemake's modified version
	- Can be used as a module or stand-alone script (anticipate more of the latter in future)
- Cutadapt parameters now directly controlled via the config
	- Default nextseq-trimming
	- Default adapter trimming
- Cluster configs converted to yaml
- Custom regex_capture for fastqs is now possible
- Conda channels reordered (faster build times)

1.4.2 (04/02/2020)
------------------
- Bugfix for regex capture w/ new file naming convention from AGC. Allows previous and current file names e.g. `_R1.fastq.gz` and `_R1_001.fastq.gz`

1.4.1 (03/25/2020)
------------------
- Refactored align_deliverables_fastq_screen (addresses errors from interference w/ parallel rule executions)
- Bugfix when order of samples in DESeqDataSet does not match samplesheet
- Bugfix for issue #20
- Cluster resource requirement tuning

1.4.0 (03/04/2020)
------------------
- Converted to singularity
    - Profiles are set-up for snakemake to use singularity by default
    - Old conda environments still work with alternative profiles, e.g. `profile-comp5-6-conda`
    - Images are auto-built by bfxcore github repositories connected to DockerHub.
- Profile added for UM-AGC (first iteration)
- Bugfix - Cutadapt was running in single-end mode in all cases. This sometimes caused issues with RSEM/STAR when resulting fastqs were out-of-sync. Now `paired_end_mode` is a cutadapt option specified in the configuration
- Bugfix - PCA plotting fix in v1.1.0 was working inconsistently.
- RScript logging is improved, though still unable to capture errors as well as we'd like.

1.3.0 (01/03/2020)
------------------
- Documentation overhaul
    - Most documents in markdown
    - All documents up-to-date
    - Made up of README, Troubleshooting doc, and several examples covering different scenarios
- Created built-in TestData genome build for watermelon_init
    - Special case - not listed in genome_references.yaml, but still usable as genome_build argument to watermelon_init
    - Sets up human chr22 references included with the repository
    - For use with the simulated data also included in the repository
    - Sets up a config perfectly compatible with README example
- Now possible to completely delete diffex portion from config, to run alignment, QC, & feature counting only
- Added functional test running snakemake with example data
- Travis-CI set up and working for repository
- Config validation bugfix - dashes no longer allowed in phenotype values
- Bugfix for unwanted trimming of column names in combined matrices
- Bugfixes for purely numeric samplenames
- Workaround to ensure watermelon_init always uses /nfs/med-bfx... over /ccmb/BioinfCore/... paths
- For consistency, all references to 'sample_description_file' changed to 'samplesheet'
- Fastqs now managed with InputFileManager (May be renamed to InputFastqManager in the future)

1.2.0 (11/14/2019)
------------------
- Added boxplot of non-normalized values
- All automated tests now working
- Added a config parameter count_min_cutoff for diffex read-count filtering threshold
- Changed summary's annotation column to entrezgene_id, since these are what we care about
- Added multiqc_rsem.txt to deliverables - this is a table with accurate alignment statistics from a single source (rsem), since multiQC isn't being very 'smart' about what it's reporting

1.1.1 (10/31/2019)
------------------
- Bugfix - Constrain snakemake to version 5.6.0 in watermelon conda env - last version supporting workflow.overwrite_configfile

1.1.0 (10/29/2019)
------------------
- Annotation overhaul
    - Ability to use GTFs from various sources as references
    - Now makes use of mapping table (source agnostic), matching gene_id field from GTF to entrezgene ID, gene symbol, and description
    - Added Rscript (ensembl_biomaRt_mapping.R) to generate said mapping table for ENSEMBL GTFs using biomaRt
    - Genome references replaced - now use ENSEMBL data - excellent versioning / standardization
    - Old genome references made compatible with new system (except E. coli, C. elegans) - still available as genome_references_old
    - annotate.py overhauled/simplified/generalized
    - annotate.py now used for both combined counts and deseq2 results
- Error logging improvements
    - Log file attached to email if errors occurred
    - Immediate email notification when first error occurs, pipeline will continue to run
- Annotated combined count matrices (from RSEM - FPKM, TPM) added to deliverables
- Gene lists are placed in one location in deliverables instead of separate folders (annotated, excel)
- References to watermelon_init are now consistent within documentation
- Fixes to make rule names & rule file names match
- Watermelon_init no longer preserves ownership/perms when copying Watermelon to project folder
- Removed strandedness option - option was accepted but not actually used by RSEM/STAR
- Added workaround for limit of 6 symbols in PCA plots
- Added guards to plotting scripts for when nrow(counts) < desired_top_n
- Added check for sample column during config validation
- Added watermelon_version to config and warning if config val != currently running val
- Added genome_build validation during watermelon_init - replaces hardcoded options
- Added ability to skip config validation
- Groundwork laid for running watermelon on GreatLakes cluster (slurm)

1.0.1 (10/22/2019)
------------------
- Bugfix - one variable wasn't renamed alongside others, causing contrast params to evaluate NULL

1.0.0 (08/22/2019)
------------------
- Watermelon (seedless)
- Replaced HISAT2/Stringtie with RSEM/STAR
- Removed ballgown
- Converted rules to use conda instead of modules
- Moved functionality from watermelon bash script directly into snakefile (eliminated watermelon.sh)
- Added the use of a snakemake profile for running on comp5/6
- Implemented samplesheet CSV to be used alongside config file
- Config validation includes samplesheet validation (modified to work with CSV input) and schema-based validation of config file
- DESeq2 parameters are directly listed in the diffex portion of the config file, offering more control over how these are run
- DESeq2 monolithic script separated into counts, init, contrasts
    - counts is run once per pipeline invocation (for a given set of alignment outputs)
    - init is run once per model, contrasts depend on this
    - contrasts is run for each contrast for a given model, i.e. they all use the same DESeq2Dataset
- Utilized snakemake's script directive, enabling snakemake S4 object to be passed directly to RScripts
- Pinned specific versions in the rule-specific conda envs, added this output to run_info deliverable
- Added output of count matrices (all samples - counts, TPM, FPKM)
- Enabled/repaired skipping of read trimming if trimming_options not set in config


0.3.6x (MM/DD/YYYY)
------------------
- Replaced Tophat2 with HISAT2; removed bbmap.
- Replaced HTSeq with stringtie; consequent renumbering of outputs
  - Added new required config value alignment_option: read_length (and set
    default to 50)
- Replaced Tuxedo/CuffDiff with Ballgown.
- Added a stand-alone Snakefile, hisat2_index.smk, to generate HISAT2 indices as necessary
- Added a top-level conda environment for watermelon
  - Upgraded Python 3.6.6, Snakemake 5.3.0, pandas (0.23.4)
- Upgraded MultiQC to 1.6 (and adjusted to use conda environment)
- Adjusted config:library_type to accept
  - fr-unstranded
  - unstranded
  - forward_reverse
  - fr-firststrand
  - reverse_forward
  - fr-secondstrand
- Throttled fastqc to avoid Java memory overallocation
- Adjusted watermelon to enable "in-flight" dry-run/dag (executed in the
  directory of a job currently in-progress).
- Added dm6 support
- Removed checksum logic in anticipation of improved sample description/
  comparison model
- Refactor diffex plots
  - DESeq2 and ballgown use the same plotting script, which requires an RData
    object from the respective diffex scripts.
  - Move plotting out of deseq2_diffex.R script
- A diffex.yaml conda environment contains all libraries needed for DESeq2, ballgown,
  and plotting. It is thus used by both diffex scripts as well as for plotting.



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
