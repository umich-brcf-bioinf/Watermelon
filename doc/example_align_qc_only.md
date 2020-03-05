# Example - Alignments, feature counts, & QC Only

Sometimes differential expression analysis is not desired, but alignment, feature counting, & QC are. This example will show how to achieve that

## The samplesheet

For this scenario, sample information relevant for differential expression is not needed, so it's possible to have a samplesheet with only the required `sample` column, like so:

    contents of samplesheet.csv
    ------------------------
    sample
    Sample_222221
    Sample_222222
    Sample_222223
    Sample_222224
    Sample_222225
    Sample_222226
    Sample_222227
    Sample_222228
    Sample_222229
    Sample_222230

Watermelon_init will still run just fine:

    #running watermelon_init
    watermelon_init --genome_build GRCh38 --job_suffix _20190823 --sample_sheet samplesheet.csv /path/to/Run_2248

Then instead of adjusting the diffex portion of the config, we can delete it entirely. This way, everything except the differential expression portions will run. Calling snakemake can be done in the same manner as the other examples:

    # Singularity must be available to snakemake, for environment management under the hood
    module load singularity/3.5.2
    snakemake --configfile config_20190823.yaml --snakefile Watermelon/rnaseq.snakefile --profile Watermelon/config/profile-greatlakes

As an additional option, it is possible to achieve the same effect leaving the config file as-is, and nullifying the diffex section directly from the command line, by adding `--config diffex=None` to the invocation like so:

    # Singularity must be available to snakemake, for environment management under the hood
    module load singularity/3.5.2
    snakemake --configfile config_20190823.yaml --snakefile Watermelon/rnaseq.snakefile --profile Watermelon/config/profile-greatlakes --config diffex=None
