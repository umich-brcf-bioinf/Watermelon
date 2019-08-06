rule align_rsem_star_combined_count_matrices:
    input:
        genes = expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results',
                sample=config['samples']),
        isoforms = expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results',
                      sample=config['samples'])
    output:
        gene = ALIGNMENT_DIR + '05-combine_counts/gene_{metric}.txt',
        isoform = ALIGNMENT_DIR + '05-combine_counts/isoform_{metric}.txt'
    log:
        ALIGNMENT_DIR + '05-combine_counts/.log/combine_counts_{metric}.log'
    params:
        genes_input_path = ALIGNMENT_DIR + '04-rsem_star_align/*.genes.results',
        isoforms_input_path = ALIGNMENT_DIR + '04-rsem_star_align/*.isoforms.results',
        metric = '{metric}'
    shell:'''(
python {WATERMELON_SCRIPTS_DIR}combine.py --input_path '{params.genes_input_path}' \
--output_file {output.gene} \
-c {params.metric} \
--id_columns gene_id

python {WATERMELON_SCRIPTS_DIR}combine.py --input_path '{params.isoforms_input_path}' \
--output_file {output.isoform} \
-c {params.metric} \
--id_columns transcript_id,gene_id

)2>&1 | tee {log}
    '''
