rule tuxedo_annotate:
    input:
        genome_checksum = CONFIG_CHECKSUMS_DIR + "config-genome.watermelon.md5",
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        gene_diff_exp = TUXEDO_DIR + "03-flag/{pheno}/{pheno}_gene.flagged.txt",
        isoform_diff_exp = TUXEDO_DIR + "03-flag/{pheno}/{pheno}_isoform.flagged.txt",
        entrez_gene_info = "references/entrez_gene_info"
    output:
        gene_annot = TUXEDO_DIR + "04-annotate/{pheno}/{pheno}_gene.flagged.annot.txt",
        isoform_annot = TUXEDO_DIR + "04-annotate/{pheno}/{pheno}_isoform.flagged.annot.txt"
    params:
        output_dir = TUXEDO_DIR + "04-annotate/{pheno}",
        genome = config["genome"]
    log:
        TUXEDO_DIR + "04-annotate/.log/{pheno}_annotate.log"
    shell:
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_annotate.py "
        " -i {input.entrez_gene_info} "
        " -e {input.gene_diff_exp} "
        " -g {params.genome} "
        " -o {params.output_dir} "
        " 2>&1 | tee {log} && "

        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_annotate.py "
        " -i {input.entrez_gene_info} "
        " -e {input.isoform_diff_exp} "
        " -g {params.genome} "
        " -o {params.output_dir} "
        " 2>&1 | tee >>{log} "
