rule align_combine_counts_to_matrices:
    input:
        genes = expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results',
            sample=SAMPLESHEET.index),
        isoforms = expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results',
            sample=SAMPLESHEET.index)
    output:
        gene = ALIGNMENT_DIR + '05-combine_counts/gene_{metric}.txt',
        isoform = ALIGNMENT_DIR + '05-combine_counts/isoform_{metric}.txt'
    log:
        JOB_LOG_DIR + 'align_combine_counts_to_matrices_{metric}.log'
    params:
        project_name = config['report_info']['project_name'],
        genes_input_path = ALIGNMENT_DIR + '04-rsem_star_align/*.genes.results',
        isoforms_input_path = ALIGNMENT_DIR + '04-rsem_star_align/*.isoforms.results',
        metric = '{metric}'
    resources: mem_mb=4000
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
