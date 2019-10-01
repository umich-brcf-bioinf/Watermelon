rule deseq2_annotation:
   input:
       diffex_file= DIFFEX_DIR + "deseq2/gene_lists/{model_name}/{contrast}.txt",
       gene_info = config['references']['annotation_tsv']
   output:
       DIFFEX_DIR + "deseq2/annotated/{model_name}/{contrast}.annot.txt"
   params:
       input_idx = 'id',
       mapping_idx = 'gene_id'
   shell:
       "python {WATERMELON_SCRIPTS_DIR}/annotate.py "
            "-m {input.gene_info} "
            "-i {input.diffex_file} "
            "-o {output} "
            "--mapping_idx {params.mapping_idx} "
            "--input_idx {params.input_idx}"
