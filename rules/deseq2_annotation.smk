rule deseq2_annotation:
   input:
       diffex_file= DIFFEX_DIR + "deseq2/gene_lists/{factor_name}/{contrast}.txt",
       gene_info = config['references']['entrez_gene_info']
   output:
       DIFFEX_DIR + "deseq2/annotated/{factor_name}/{contrast}.annot.txt"
   params:
       genome = config["genome"]
   shell:
       "python {WATERMELON_SCRIPTS_DIR}/deseq2_annotate.py "
       "    -i {input.gene_info} "
       "    -e {input.diffex_file} "
       "    -g {params.genome} "
       "    -o {output}"
