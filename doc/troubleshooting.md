# Troubleshooting

## Locked directory

A locked directory happens sometimes when snakemake is quit in the middle of a job. The error may say something like:

    Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory:
    /ccmb/BioinfCore/SoftwareDev/projects/Watermelon_spike/Tronson_RS1_tests/test_data/outputs/012017/test_run_watermelon/Demo_test1/test/analysis_fewer_replicates
    If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the --unlock argument.

As noted, this can be resolved by appending `--unlock` to your previous attempted snakemake command. After unlocking is successfully completed, then you can run the original command again and it should work.

## Rerun-incomplete

Another error that sometimes occurs if a job is ended unexpectedly, is the rerun-incomplete error. The solution to this depends on state of the files in question. The message may say something like:

    The files below seem to be incomplete. If you are sure that certain files are not incomplete, mark them as complete with

      snakemake --cleanup-metadata <filenames>`

    To re-generate the files rerun your command with the --rerun-incomplete flag.
    Incomplete files:
    alignment_results/03-fastqc_reads/reads_fastq.done
    alignment_results/03-fastqc_reads/Sample_61492_trimmed_R1_fastqc.html

As described in the error message, if you're sure that the files are complete, then it's possible to mark them as complete by adding `--cleanup-metadata` to the command line options. Alternatively, if it's suspected that the files are truly incomplete, e.g. if the job might have been killed prematurely, then you can similarly append `--rerun-incomplete` to the snakemake command, and it will delete those files and re-run their jobs.
