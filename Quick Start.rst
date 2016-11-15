
-----------
Quick Start
-----------

1. **Add 'bfx_modules' to your bashrc:**

Modify .bashrc to add this line: 
::
$ module use /nfs/med-bfx-common/software/bfx_modules
  
2. **Download the example data:**
::

$ git clone <umich-brcf-bioinf/Watermelon.git...? sample_data

3. **Run 'watermelon-init' to setup the analysis**
  Requires: 
* genome build
* path to sample directories with multiplexed reads
* project tag
::

  $ watermelon-init mm10 /ccmb/BioinfCore/DNASeqCore/Run_1639/ksinger/Run_1639/ksinger Singer_RS1_DietCell

4. **Review 'watermelon-init' output:**
  Generates three directories: 
* inputs : Symlinks to the samples dirs of multiplexed sequences
* analysis-project_tag  : config.yaml file (which needs to be set-up to run the analysis)
* deliverables-project_tag :Results of the analysis
::

  $ ls
    inputs/
      Sample_61483/
      Sample_61484/
      Sample_61490/
      Sample_61491/
      Sample_61501/
      Sample_61502/
      Sample_61505/
      Sample_61506/
    
    analysis-Mouse_RS1_Condition1/
        Mouse_RS1_Condition1_config.yaml
        references/
        
    deliverables-Mouse_RS1_Condition1/

5. **Setup Watermelon analysis**
::

  $ cd analysis-Mouse_RS1_Condition1

* Edit the config file (Mouse_RS1_Condition1_config.yaml). 
* Add comparison details; set trimming, alignment options, and fold change threshold.

6. **Run Watermelon**
::
  $ watermelon -c Mouse_RS1_Condition1_config.yaml

7. **Reviewing output files**
  Watermelon generates the following result files:
  
* fastqc reports of raw reads
* fastqc reports of aligned bams 
* alignment summary stats

* cuffdiff report excel files
* cummerbund plots
* tab-delimited file of raw counts of fragments per sample


::

  $ ls deliverables-Mouse_RS1_Condition1/
  
      diffex/
        cummerbund_plots/
        LVM_v_LVF.xlsx  
        VVM_v_LLF.xlsx
        LLF_v_LVF_v_LVM_v_VVM_repRawCounts.txt
        
      qc/
      raw_reads_fastqc/
      aligned_reads_fastqc/
      align_summary.txt

**8. How to restart Watermelon?**

If you quit Watermelon in the middle of a run, and want to restart it, Watermelon will automatically restart from the point you exited. If partially created files exist, Watermelon will prompt you to use the following commands to clean-up partially created directories and files and restart analysis. 
::

  --unlock              Remove a lock on the working directory.
  --cleanup-metadata FILE [FILE ...], --cm  Cleanup the metadata of given files. 


If you made changes to the config.yaml (e.g. Added new samples, comparisons, changes to trimming options, etc.), then Watermelon will restart the workflow in a way that meets the workflow’s logical requirements. 

For example, if you added two new samples to the config, Watermelon, will run the Concat, QC and alignment steps for these two samples, and then re-run the entire differential expression analysis. If you just added a new comparison (for existing samples), then it will simply re-run the differential expression analysis steps.

