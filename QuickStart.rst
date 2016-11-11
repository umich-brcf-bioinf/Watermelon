-----------
Quick Start
-----------

1. **Add 'bfx_modules' to your bashrc**

Modify .bashrc to add this line: 
::

  $ module use /nfs/med-bfx-common/software/bfx_modules
  
2. **Run Watermelon-init**
watermelon-init sets up your anlaysis. Provide the following: 
  * Path to sample directories with multiplexed sequence reads (e.g./ccmb/BioinfCore/DNASeqCore/Run_1639/ksinger/Run_1639/ksinger)
  * genome build (e.g. mm10)
  * project tag (e.g. Singer_RS1_DietCell)
::

  $ watermelon-init mm10 /ccmb/BioinfCore/DNASeqCore/Run_1639/ksinger/Run_1639/ksinger Singer_RS1_DietCell

This will generate three directories: 
    * inputs
    * analysis-Singer_RS1_DietCell
    * deliverables-Singer_RS1_DietCell



  
