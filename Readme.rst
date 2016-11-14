
**Watermelon directory structure**

Input directory
* Contains the multiplexed raw reads
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
* Contains config file, and sub-directories for each step of the analysis
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
* Contains the analysis results
::

  deliverables-Mouse_RS1_Condition1
      diffex/
      qc/
