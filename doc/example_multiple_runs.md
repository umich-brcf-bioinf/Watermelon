# Example - Running with multiple sequencing runs from UMich Advanced Genomics Core

The walkthrough in the readme provided an example that should work out-of-the-box, where things were simplified for illustrative purposes. The current example goes a bit further, highlighting a few more details which might be useful when running with multiple sequencing runs.

## Before starting:

Just as in the readme walkthrough, you'll run watermelon_init to set up the project directory. Also like before, the watermelon conda environment must be activated and the path to Watermelon should be in your PATH environment variable in order to run watermelon_init.

## Multiple sequencing runs, different biological samples

There are a couple of scenarios worth considering. The first is the simplest - when you want to combine sequencing runs in a workflow, but the samples included in each run are distinct. In this case, you can provide watermelon_init with multiple input_run_dirs, and each sample will be mapped to its appropriate input directory in the samplesheet.

    watermelon_init.py --genome_build GRCh38 --project_id 20190822 --type align_qc --input_run_dirs /path/to/Run_1234 /path/to/Run_1256

## Same biological samples over multiple sequencing runs

The next scenario is a bit more involved in setting up. Here we'll have a subset of biological samples that were sequenced across two separate runs. In this case, we want to merge the samples across the runs, so that we can effectively achieve greater sequencing depth by concatenating the samples' fastq files before alignment. I'll illustrate how to handle this scenario.

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


To achieve the result that we want, we should manually move the desired files from each sample into a single directory, so that these fastq files will be concatenated before alignment.

Here's some code that would work for the above example. Note that it uses the perl-based rename program, not the util-linux one. To create & activate a perl_rename conda environment:

    conda create -n perl_rename -c bioconda rename
    conda activate perl_rename

With rsync and perl_rename, we can use the following to merge our input run directories. Note that with rsync, we use the '-b' flag to prevent overwriting and a '--suffix' to control the names. This isn't an issue with this example, but this is something to be cognizant of. Also note the trailing slash on the source directory `run_1256/`, it's important to match the syntax shown below for the correct behavior.

    # First merge the directories with rsync.
    rsync -ab --suffix _run_1256 run_1256/ run_1234
    # Then do dry-run of rename
    find run_1234/ -type f -name "*.fastq.gz_run_1256" -exec rename -ne 's#(.*)/(.*)_run_1256$#$1/run_1256_$2#' {} \;
    # Finally the rename
    find run_1234/ -type f -name "*.fastq.gz_run_1256" -exec rename -e 's#(.*)/(.*)_run_1256$#$1/run_1256_$2#' {} \;

After this, we would have the following output directory structure:

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
    │   ├── run_1256_Sample_358629_L001_R1.fastq.gz
    │   ├── run_1256_Sample_358629_L001_R2.fastq.gz
    │   ├── run_1256_Sample_358629_L002_R1.fastq.gz
    │   ├── run_1256_Sample_358629_L002_R2.fastq.gz
    │   ├── Sample_358629_L003_R1.fastq.gz
    │   ├── Sample_358629_L003_R2.fastq.gz
    │   ├── Sample_358629_L004_R1.fastq.gz
    │   └── Sample_358629_L004_R2.fastq.gz
    ├── Sample_358630
    │   ├── run_1256_Sample_358630_L001_R1.fastq.gz
    │   ├── run_1256_Sample_358630_L001_R2.fastq.gz
    │   ├── run_1256_Sample_358630_L002_R1.fastq.gz
    │   ├── run_1256_Sample_358630_L002_R2.fastq.gz
    │   ├── Sample_358630_L003_R1.fastq.gz
    │   ├── Sample_358630_L003_R2.fastq.gz
    │   ├── Sample_358630_L004_R1.fastq.gz
    │   └── Sample_358630_L004_R2.fastq.gz
    ├── Sample_358631
    │   ├── run_1256_Sample_358631_L001_R1.fastq.gz
    │   ├── run_1256_Sample_358631_L001_R2.fastq.gz
    │   ├── run_1256_Sample_358631_L002_R1.fastq.gz
    │   ├── run_1256_Sample_358631_L002_R2.fastq.gz
    │   ├── Sample_358631_L003_R1.fastq.gz
    │   ├── Sample_358631_L003_R2.fastq.gz
    │   ├── Sample_358631_L004_R1.fastq.gz
    │   └── Sample_358631_L004_R2.fastq.gz
    ├── Sample_358632
    │   ├── run_1256_Sample_358632_L001_R1.fastq.gz
    │   ├── run_1256_Sample_358632_L001_R2.fastq.gz
    │   ├── run_1256_Sample_358632_L002_R1.fastq.gz
    │   ├── run_1256_Sample_358632_L002_R2.fastq.gz
    │   ├── Sample_358632_L003_R1.fastq.gz
    │   ├── Sample_358632_L003_R2.fastq.gz
    │   ├── Sample_358632_L004_R1.fastq.gz
    │   └── Sample_358632_L004_R2.fastq.gz
    ├── Sample_358633
    │   ├── run_1256_Sample_358633_L001_R1.fastq.gz
    │   ├── run_1256_Sample_358633_L001_R2.fastq.gz
    │   ├── run_1256_Sample_358633_L002_R1.fastq.gz
    │   ├── run_1256_Sample_358633_L002_R2.fastq.gz
    │   ├── Sample_358633_L003_R1.fastq.gz
    │   ├── Sample_358633_L003_R2.fastq.gz
    │   ├── Sample_358633_L004_R1.fastq.gz
    │   └── Sample_358633_L004_R2.fastq.gz
    └── Sample_358634
        ├── run_1256_Sample_358634_L001_R1.fastq.gz
        ├── run_1256_Sample_358634_L001_R2.fastq.gz
        ├── run_1256_Sample_358634_L002_R1.fastq.gz
        ├── run_1256_Sample_358634_L002_R2.fastq.gz
        ├── Sample_358634_L003_R1.fastq.gz
        ├── Sample_358634_L003_R2.fastq.gz
        ├── Sample_358634_L004_R1.fastq.gz
        └── Sample_358634_L004_R2.fastq.gz

## Watermelon_init
After this, we can run watermelon_init, supplying it with the merged input run directory.

    #In the project directory, run watermelon_init
    watermelon_init.py --genome_build GRCh38 --project_id 20190822 --type align_qc --input_run_dirs /path/to/Run_1234


## Running the pipeline

After looking over the config, running the pipeline is the same as in the main example:

    # Start a screen session (for persistence over ssh):
    screen -S watermelon_20190822
    # Activate the conda environment:
    conda activate watermelon
    # Dry-run to validate the config and check the execution plan:
    snakemake --dryrun --printshellcmds --configfile config_20190822.yaml --snakefile Watermelon/align_qc.smk

Still in the project directory, now ready to run the pipeline - run on bfx-comp5/6 (notice the profile):

    # Singularity must be available to snakemake, for environment management under the hood
    module load singularity
    snakemake --configfile config_20190822.yaml --snakefile Watermelon/align_qc.smk --profile Watermelon/config/profile-comp5-6
