-----------
Quick Start
-----------

1. **Add 'bfx_modules' to your bashrc**

  Modify .bashrc to add this line: 
::
  $ module use /nfs/med-bfx-common/software/bfx_modules
  
2. **Run watermelon-init to setup the analysis**
  Requires: 
    * genome build
    * path to sample directories with multiplexed reads
    * project tag
::

  $ watermelon-init mm10 /ccmb/BioinfCore/DNASeqCore/Run_1639/ksinger/Run_1639/ksinger Singer_RS1_DietCell

3. **Reviewing the watermelon-init output**
  This command will generate three directories: 
    * inputs : contains symlinks to the samples dirs of multiplexed sequences
    * analysis-project_tag  : This is where the analysis will be run. Contains the config.yaml file which needs to be set-up to run the analysis.
    * deliverables-project_tag (e.g. ) :Contains the results of the analysis

  * Input directory
  
::
  inputs/
    Sample_61483/
    Sample_61484/
    Sample_61490/
    Sample_61491/
    Sample_61501/
    Sample_61502/
    Sample_61505/
    Sample_61506/
 
  Analysis directory
::
  analysis-Mouse_RS1_Condition1/
  
      Mouse_RS1_Condition1_config.yaml
      
      01-raw_reads/
      02-cutadapt/
      03-fastqc_reads/
      04-tophat/
      05-fastqc_align/
      06-qc_metrics/
      07-htseq/
      08-cuffdiff/
      09-diffex_flip/
      10-diffex_flag/
      11-annotate_diffex_flag/
      12-group_replicates/
      13-cummerbund/
      14-diffex_split/
      15-diffex_excel/
      16-deliverables/
      config_checksums/
      logs/
      references/
      
  Deliverables directory
::
  deliverables-Mouse_RS1_Condition1
      diffex/
      qc/


3 **Set up watermelon analysis**
  * cd into 'analysis' directory. 
  * Edit the config file (Mouse_RS1_Condition1_config.yaml). Add comparison details; set trimming, alignment options, and fold change threshold.

4. **Run watermelon**
::
  $ watermelon
