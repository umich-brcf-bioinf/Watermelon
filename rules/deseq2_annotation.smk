rule deseq2_annotation:
   input:
       diffex_file= DESEQ2_DIR + "03-deseq2_diffex/gene_lists/{phenotype}/{comparison}.txt",
       gene_info = "references/entrez_gene_info",
   output:
       DESEQ2_DIR + "04-annotation/{phenotype}/{comparison}.annot.txt",
   params:
       genome = config["genome"],
       output_dir = DESEQ2_DIR + "04-annotation/{comparison}",
   shell:
       "rm -rf {params.output_dir}.tmp && "
       "mkdir -p {params.output_dir}.tmp && "
       "python {WATERMELON_SCRIPTS_DIR}/deseq2_annotate.py "
       "    -i {input.gene_info} "
       "    -e {input.diffex_file} "
       "    -g {params.genome} "
       "    -o {output}.tmp && "
       "mv {output}.tmp {output} "
