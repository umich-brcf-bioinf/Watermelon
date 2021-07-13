rule deseq2_annotation:
   input:
       diffex_file = DIFFEX_DIR + "{model_name}/gene_lists/{contrast}.txt",
       gene_info = config['references']['annotation_tsv']
   output:
       annot_results = DIFFEX_DIR + "{model_name}/annotated/{contrast}.annot.txt"
   params:
        project_name = config['report_info']['project_name'],
        input_idx = 'gene_id',
        mapping_idx = 'gene_id'
   shell:
       "python {WATERMELON_SCRIPTS_DIR}/annotate.py "
            "-m {input.gene_info} "
            "-i {input.diffex_file} "
            "-o {output} "
            "--mapping_idx {params.mapping_idx} "
            "--input_idx {params.input_idx}"
