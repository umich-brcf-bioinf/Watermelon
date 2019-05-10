rule ballgown_annotation:
    input:
        genome_checksum = CONFIG_CHECKSUMS_DIR + "config-genome.watermelon.md5",
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        gene_diff_exp = BALLGOWN_DIR + '01-ballgown_diffex/gene_lists/{phenotype}/{comparison}_gene.txt',
        isoform_diff_exp = BALLGOWN_DIR + '01-ballgown_diffex/gene_lists/{phenotype}/{comparison}_isoform.txt',
        entrez_gene_info = "references/entrez_gene_info"
    output:
        gene_annot = BALLGOWN_DIR + "02-annotate/{phenotype}/{comparison}_gene.annot.txt",
        isoform_annot = BALLGOWN_DIR + "02-annotate/{phenotype}/{comparison}_isoform.annot.txt"
    params:
        output_dir = BALLGOWN_DIR + "02-annotate/{phenotype}",
        genome = config["genome"]
    log:
        BALLGOWN_DIR + "02-annotate/.log/{phenotype}_{comparison}_annotate.log"
    shell:
        '''
        python {WATERMELON_SCRIPTS_DIR}/ballgown_annotate.py \
         -i {input.entrez_gene_info} \
         -e {input.gene_diff_exp} \
         -g {params.genome} \
         -o {params.output_dir} \
         2>&1 | tee {log}

        python {WATERMELON_SCRIPTS_DIR}/ballgown_annotate.py \
         -i {input.entrez_gene_info} \
         -e {input.isoform_diff_exp} \
         -g {params.genome} \
         -o {params.output_dir} \
         2>&1 | tee {log}
        '''
