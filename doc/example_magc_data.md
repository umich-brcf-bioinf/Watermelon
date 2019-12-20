# Example - running with data from MAGC

The walkthrough in the readme provided an example that should work out-of-the-box, where things were simplified for illustrative purposes. There are a few more details to pay attention to when running with real data. These will be highlighted in this example. One major

## Before starting:

Just as in the readme walkthrough, you'll run watermelon_init to set up the project directory. Also like before, the watermelon conda environment must be activated and the path to Watermelon/bin should be in your PATH environment variable in order to run watermelon_init.

## Fastq directories for watermelon_init

Say for this example that we have a subset of biological samples that were sequenced across two separate runs. I'll highlight again here the requirements for the input structure, and also illustrate the capability of watermelon_init to handle this scenario.

The input directories given as command-line arguments to watermelon_init should themselves contain subdirectories named after the samples. In this example we'll have two such top-level directories, one for each sequencing run.

Run_1234 contains 4 lanes of paired-end sequencing data for half of the samples, and 2 lanes of paired-end data for the second half.

Run 1256 contains an additional 2 lanes of paired-end sequencing data for the second half of samples.

    contents/structure of Run_1234/
    ---------------------
    Run_1234/
    ├── Sample_358623
    │   ├── Sample_358623_L001_R1.fastq.gz
    │   ├── Sample_358623_L001_R2.fastq.gz
    │   ├── Sample_358623_L002_R1.fastq.gz
    │   ├── Sample_358623_L002_R2.fastq.gz
    │   ├── Sample_358623_L003_R1.fastq.gz
    │   ├── Sample_358623_L003_R2.fastq.gz
    │   ├── Sample_358623_L004_R1.fastq.gz
    │   └── Sample_358623_L004_R2.fastq.gz
    ├── Sample_358624
    │   ├── Sample_358624_L001_R1.fastq.gz
    │   ├── Sample_358624_L001_R2.fastq.gz
    │   ├── Sample_358624_L002_R1.fastq.gz
    │   ├── Sample_358624_L002_R2.fastq.gz
    │   ├── Sample_358624_L003_R1.fastq.gz
    │   ├── Sample_358624_L003_R2.fastq.gz
    │   ├── Sample_358624_L004_R1.fastq.gz
    │   └── Sample_358624_L004_R2.fastq.gz
    ├── Sample_358625
    │   ├── Sample_358625_L001_R1.fastq.gz
    │   ├── Sample_358625_L001_R2.fastq.gz
    │   ├── Sample_358625_L002_R1.fastq.gz
    │   ├── Sample_358625_L002_R2.fastq.gz
    │   ├── Sample_358625_L003_R1.fastq.gz
    │   ├── Sample_358625_L003_R2.fastq.gz
    │   ├── Sample_358625_L004_R1.fastq.gz
    │   └── Sample_358625_L004_R2.fastq.gz
    ├── Sample_358626
    │   ├── Sample_358626_L001_R1.fastq.gz
    │   ├── Sample_358626_L001_R2.fastq.gz
    │   ├── Sample_358626_L002_R1.fastq.gz
    │   ├── Sample_358626_L002_R2.fastq.gz
    │   ├── Sample_358626_L003_R1.fastq.gz
    │   ├── Sample_358626_L003_R2.fastq.gz
    │   ├── Sample_358626_L004_R1.fastq.gz
    │   └── Sample_358626_L004_R2.fastq.gz
    ├── Sample_358627
    │   ├── Sample_358627_L001_R1.fastq.gz
    │   ├── Sample_358627_L001_R2.fastq.gz
    │   ├── Sample_358627_L002_R1.fastq.gz
    │   ├── Sample_358627_L002_R2.fastq.gz
    │   ├── Sample_358627_L003_R1.fastq.gz
    │   ├── Sample_358627_L003_R2.fastq.gz
    │   ├── Sample_358627_L004_R1.fastq.gz
    │   └── Sample_358627_L004_R2.fastq.gz
    ├── Sample_358628
    │   ├── Sample_358628_L001_R1.fastq.gz
    │   ├── Sample_358628_L001_R2.fastq.gz
    │   ├── Sample_358628_L002_R1.fastq.gz
    │   ├── Sample_358628_L002_R2.fastq.gz
    │   ├── Sample_358628_L003_R1.fastq.gz
    │   ├── Sample_358628_L003_R2.fastq.gz
    │   ├── Sample_358628_L004_R1.fastq.gz
    │   └── Sample_358628_L004_R2.fastq.gz
    ├── Sample_358629
    │   ├── Sample_358629_L003_R1.fastq.gz
    │   ├── Sample_358629_L003_R2.fastq.gz
    │   ├── Sample_358629_L004_R1.fastq.gz
    │   └── Sample_358629_L004_R2.fastq.gz
    ├── Sample_358630
    │   ├── Sample_358630_L003_R1.fastq.gz
    │   ├── Sample_358630_L003_R2.fastq.gz
    │   ├── Sample_358630_L004_R1.fastq.gz
    │   └── Sample_358630_L004_R2.fastq.gz
    ├── Sample_358631
    │   ├── Sample_358631_L003_R1.fastq.gz
    │   ├── Sample_358631_L003_R2.fastq.gz
    │   ├── Sample_358631_L004_R1.fastq.gz
    │   └── Sample_358631_L004_R2.fastq.gz
    ├── Sample_358632
    │   ├── Sample_358632_L003_R1.fastq.gz
    │   ├── Sample_358632_L003_R2.fastq.gz
    │   ├── Sample_358632_L004_R1.fastq.gz
    │   └── Sample_358632_L004_R2.fastq.gz
    ├── Sample_358633
    │   ├── Sample_358633_L003_R1.fastq.gz
    │   ├── Sample_358633_L003_R2.fastq.gz
    │   ├── Sample_358633_L004_R1.fastq.gz
    │   └── Sample_358633_L004_R2.fastq.gz
    └── Sample_358634
        ├── Sample_358634_L003_R1.fastq.gz
        ├── Sample_358634_L003_R2.fastq.gz
        ├── Sample_358634_L004_R1.fastq.gz
        └── Sample_358634_L004_R2.fastq.gz

    12 directories, 72 files


    contents/structure of Run_1256
    ------------------------------
    Run_1256/
    ├── Sample_358629
    │   ├── Sample_358629_L001_R1.fastq.gz
    │   ├── Sample_358629_L001_R2.fastq.gz
    │   ├── Sample_358629_L002_R1.fastq.gz
    │   └── Sample_358629_L002_R2.fastq.gz
    ├── Sample_358630
    │   ├── Sample_358630_L001_R1.fastq.gz
    │   ├── Sample_358630_L001_R2.fastq.gz
    │   ├── Sample_358630_L002_R1.fastq.gz
    │   └── Sample_358630_L002_R2.fastq.gz
    ├── Sample_358631
    │   ├── Sample_358631_L001_R1.fastq.gz
    │   ├── Sample_358631_L001_R2.fastq.gz
    │   ├── Sample_358631_L002_R1.fastq.gz
    │   └── Sample_358631_L002_R2.fastq.gz
    ├── Sample_358632
    │   ├── Sample_358632_L001_R1.fastq.gz
    │   ├── Sample_358632_L001_R2.fastq.gz
    │   ├── Sample_358632_L002_R1.fastq.gz
    │   └── Sample_358632_L002_R2.fastq.gz
    ├── Sample_358633
    │   ├── Sample_358633_L001_R1.fastq.gz
    │   ├── Sample_358633_L001_R2.fastq.gz
    │   ├── Sample_358633_L002_R1.fastq.gz
    │   └── Sample_358633_L002_R2.fastq.gz
    └── Sample_358634
        ├── Sample_358634_L001_R1.fastq.gz
        ├── Sample_358634_L001_R2.fastq.gz
        ├── Sample_358634_L002_R1.fastq.gz
        └── Sample_358634_L002_R2.fastq.gz

    6 directories, 24 files

For the samples which were sequenced across multiple sequencing runs, watermelon_init will effectively merge these directories. We'll supply the path to both Run_1234 and Run_1256 as arguments to watermelon_init.

## The samplesheet

The samplesheet contains a column `sample` containing sample IDs, and additional columns with sample information which will be used for the differential expression analysis. These columns will come into play after watermelon_init is run, during adjustment of the config.

    contents of samplesheet.csv
    ------------------------
    sample,phenotype
    Sample_358623,mutA
    Sample_358624,mutA
    Sample_358625,mutA
    Sample_358626,mutA
    Sample_358627,mutB
    Sample_358628,mutB
    Sample_358629,mutB
    Sample_358630,mutB
    Sample_358631,WT
    Sample_358632,WT
    Sample_358633,WT
    Sample_358634,WT

## Running watermelon_init

    #In the project directory (where the analysis will be run and results will be stored), run watermelon_init
    watermelon_init --genome_build GRCh38 --job_suffix _20190822 --sample_sheet samplesheet.csv /path/to/Run_1234 /path/to/Run_1256

## Adjusting the config

The configuration output by watermelon_init provides a template which should be adjusted for the analysis at hand. The file `watermelon.README` also generated by watermelon_init details some of these necessary adjustments as well, so that is useful to review. For now we'll assume that only the differential expression analysis portion needs to be adjusted. Here's that portion excerpted from the config file (I removed comments to make it more readable here):

    diffex:
        linear_fold_change: 1.5
        adjustedPValue: 0.05
        count_min_cutoff: 1
        model_treatment:
            contrasts:
            -   drug_v_control
            DESeq2:
                design: ~ treatment
                factor_name: treatment
                DESeq2:
                    betaPrior: true
                results:
                    contrast: list(paste0(factorName, testName), paste0(factorName, referenceName))

To adjust this config and get it working, we'll change the name of the model, change the list of contrasts, the design, and the factor name to match our experiment. The adjusted config will have the following values for this analysis:

    diffex:
        linear_fold_change: 1.5
        adjustedPValue: 0.05
        count_min_cutoff: 1
        model_mutation:
            contrasts:
            -   mutA_v_WT
            -   mutB_v_WT
            DESeq2:
                design: ~ mutation
                factor_name: mutation
                DESeq2:
                    betaPrior: true
                results:
                    contrast: list(paste0(factorName, testName), paste0(factorName, referenceName))

## Running the pipeline

After adjusting the config, running the pipeline is the same as in the main example:

    # Start a screen session (for persistence over ssh):
    screen -S watermelon_20190822
    # Activate the conda environment:
    conda activate watermelon
    # Dry-run to validate the config and check the execution plan:
    snakemake --dryrun --printshellcmds --configfile config_20190822.yaml --snakefile Watermelon/rnaseq.snakefile

Still in the project directory, now ready to run the pipeline - run on bfx-comp5/6 (notice the profile):

    snakemake --use-conda --configfile config_20190822.yaml --snakefile Watermelon/rnaseq.snakefile --profile Watermelon/config/profile-comp5-6
