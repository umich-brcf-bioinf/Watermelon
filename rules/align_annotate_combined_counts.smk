rule align_annotate_combined_counts:
   input:
       counts_file = ALIGNMENT_DIR + "05-combine_counts/{feature}_{metric}.txt",
       gene_info = config['references']['annotation_tsv']
   output:
       ALIGNMENT_DIR + "06-annotate_combined_counts/{feature}_{metric}.annot.txt"
   params:
       input_idx = 'gene_id',
       mapping_idx = 'gene_id'
   shell:
       "python {WATERMELON_SCRIPTS_DIR}/annotate.py "
            "-m {input.gene_info} "
            "-i {input.counts_file} "
            "-o {output} "
            "--mapping_idx {params.mapping_idx} "
            "--input_idx {params.input_idx}"
