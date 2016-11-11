
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

3. **Review/edit watermelon-init output**
  This command will generate three directories: 
    * inputs : Symlinks to the samples dirs of multiplexed sequences
    * analysis-project_tag  : config.yaml file (which needs to be set-up to run the analysis)
    * deliverables-project_tag :Results of the analysis
  
::

  $ ls inputs/
  
    Sample_61483/
    Sample_61484/
    Sample_61490/
    Sample_61491/
    Sample_61501/
    Sample_61502/
    Sample_61505/
    Sample_61506/

::

  $ ls analysis-Mouse_RS1_Condition1/
  
    Mouse_RS1_Condition1_config.yaml
::

  $ ls deliverables-Mouse_RS1_Condition1/
  
      diffex/
        cummerbund_plots/
        
        LVM_v_LVF.xlsx  
        VVM_v_LLF.xlsx
        LLF_v_LVF_v_LVM_v_VVM_repRawCounts.txt
        
      qc/

3 **Set up watermelon analysis**
  * cd into 'analysis' directory. 
  * Edit the config file (Mouse_RS1_Condition1_config.yaml). Add comparison details; set trimming, alignment options, and fold change threshold.

4. **Run watermelon**
::
  $ watermelon
