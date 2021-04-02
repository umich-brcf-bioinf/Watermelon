# Example - Starting from count matrix

The initialization program `watermelon_init` has the ability to intake a count matrix file instead of a list of directories. This is useful in situations where the count matrix is already available. Similar to the other examples, the watermelon conda environment must be activated and the path to Watermelon/bin should be in your PATH environment variable in order to run watermelon_init.

    watermelon_init --genome_build GRCh38 --job_suffix _20190822 --sample_sheet samplesheet.csv --count_matrix counts_ready.txt

This will place the count matrix into the resulting configuration, and also set a few configuration details to prepare the pipeline to run from counts onward. The snakemake command to execute the pipeline is the same as in other examples - the changes within the configuration file will ensure that the count table is used to start the pipeline.

Count matrix details:
- Text file of tab-separated values
- Row for each gene, matching the annotation in genome_build
- Column for each sample with count values, with column names matching the sample sheet_file

Note:
It is important that the count matrix only contains the above information. If annotation columns are present (such as description), these must be removed before running the pipeline

    # e.g. if columns 2,3,4 have annotation information, gather all the rest
    cut -d $'\t' -f 1,5- /path/to/gene_expected_count.annot.txt > counts_ready.txt
