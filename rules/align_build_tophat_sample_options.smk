rule align_build_tophat_sample_options:
    input:
        read_stats = lambda wildcards: rnaseq_snakefile_helper.expand_read_stats_if_paired(\
                ALIGNMENT_DIR + "04-insert_size/{sample}_read_stats.txt",
                SAMPLE_READS,
                wildcards.sample)
    output:
        tophat_sample_options = ALIGNMENT_DIR + "04-tophat/{sample}/{sample}.tophat_sample_options",
    run:
        options = ''
        with open(str(output['tophat_sample_options']), 'w') as output:
            if 'read_stats' in input:
                options = rnaseq_snakefile_helper.tophat_paired_end_flags(str(input['read_stats']))
            output.write(options)
