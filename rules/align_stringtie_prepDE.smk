# ALL.append([ALIGNMENT_DIR + '06-stringtie/gene_count_matrix.csv',
#             ALIGNMENT_DIR + '06-stringtie/transcript_count_matrix.csv'])

rule align_stringtie_prepDE:
    input:
        gtfs = expand(ALIGNMENT_DIR + '06-stringtie/{sample}.gtf',
                      sample=config[SAMPLES_KEY]),
        prep_input = expand(ALIGNMENT_DIR + '06-stringtie/{sample}_prepDE_input.txt',
                            sample=config[SAMPLES_KEY]),
    output:
        prep_output = ALIGNMENT_DIR + '06-stringtie/prepDE_input.txt',
        gene_counts = ALIGNMENT_DIR + '06-stringtie/gene_count_matrix.csv',
        transcript_counts = ALIGNMENT_DIR + '06-stringtie/transcript_count_matrix.csv',
    log:
        ALIGNMENT_DIR + '06-stringtie/.log/align_stringtie_prepDE.log'
    # benchmark:
    #     'benchmarks/align_stringtie_prepDE.benchmark.txt'
    conda:
        'envs/align_hisat2_stringtie.yaml'
    params:
        length = config['alignment_options']['read_length'],
        prep_input = ' '.join(expand(ALIGNMENT_DIR + '06-stringtie/{sample}_prepDE_input.txt',
                                     sample=config[SAMPLES_KEY]))
    shell:
        '''(
        cat {params.prep_input} > {output.prep_output}

        prepDE.py \
            -i {output.prep_output} \
            -g {output.gene_counts} \
            -t {output.transcript_counts} \
            -l {params.length}
        ) 2>&1 | tee {log}'''
