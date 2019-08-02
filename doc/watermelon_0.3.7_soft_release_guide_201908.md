# Watermelon 0.3.7 soft-release instructions

Between Watermelon v0.3.6 and v0.3.7, there have been a large number of significant changes to the codebase, in addition to significant changes in how the analyst/user will launch Watermelon. The soft-release of Watermelon v0.3.7 will be used to familiarize users with these changes, and to get feedback on any points that may require attention before a full release. It is also worth mentioning that the existing documentation is outdated and a large percentage of it must be rewritten. The information gathered during the soft release will be used to determine the most important points that must be addressed in the documentation

### Getting started

First, we will need to get the develop branch of the Watermelon pipeline. I suggest creating a couple of directories - one for the code, and one for a test project

    cd /nfs/med-bfx-activeprojects/$USER
    mkdir watermelon-0.3.7-develop
    mkdir watermelon-0.3.7-testproject
    #Now get the develop version of watermelon
    cd watermelon-0.3.7-develop
    git clone --single-branch --branch develop https://github.com/umich-brcf-bioinf/Watermelon/

Since this is the first time an analyst is running watermelon seedless, it requires the preparation of a conda environment before watermelon can be run. (Hopefully anaconda/miniconda is already installed. If not, see ). The `conda env create` command below will create a conda environment in the user's anaconda installation location entitled "watermelon".
The conda environment only needs to be created at most once per release. However, the conda environment must be activated in each new terminal before the pipeline can be used (similar to loading a module).

    conda env create -f Watermelon/envs/watermelon.yaml
    #After completion, activate the conda environment
    conda activate watermelon
    #Take a look and see that your watermelon environment is activated (it will be starred)
    conda info --envs

Now go to the project directory and grab the example samplesheet.
In an actual use-case, the analyst would prepare a samplesheet according to their project plan
and place it here (in the project directory).

    cd ../watermelon-0.3.7-testproject
    cp ../watermelon-0.3.7-develop/Watermelon/config/example_sample_description.csv .

Run watermelon_init - this will gather the input fastqs, validate the samplesheet, and prepare a basic configuration file for the analysis

    ../watermelon-0.3.7-develop/Watermelon/watermelon_init-runner.py \
      --genome_build hg19 --job_suffix "_test_$(date +%F)" \
      --sample_sheet example_sample_description.csv \
      /nfs/med-bfx-dnaseqcore/Run_1562/lumeng/Run_1562/lumeng/ \
      /nfs/med-bfx-dnaseqcore/Run_1564/lumeng/Run_1564/lumeng/ \
      /nfs/med-bfx-dnaseqcore/Run_1567/lumeng/Run_1567/lumeng/
