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
    
::
  $ ls
     inputs
     analysis-Singer_RS1_DietCell
     deliverables-Singer_RS1_DietCell
    

3 **Set up watermelon analysis**

cd into analysis-project_tag directory; edit the config file: project_tag_config.yaml, 

4. **Run watermelon**
