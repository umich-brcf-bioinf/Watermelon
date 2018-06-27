rule tuxedo_flag:
    input:
        fold_change_checksum = CONFIG_CHECKSUMS_DIR + "config-fold_change.watermelon.md5",
        cuffdiff_gene_exp = TUXEDO_DIR + "02-flip/{pheno}/gene_exp.flip.diff",
        cuffdiff_isoform_exp = TUXEDO_DIR + "02-flip/{pheno}/isoform_exp.flip.diff"
    output:
        gene_flagged = TUXEDO_DIR + "03-flag/{pheno}/{pheno}_gene.flagged.txt",
        isoform_flagged = TUXEDO_DIR + "03-flag/{pheno}/{pheno}_isoform.flagged.txt",
    params:
        fold_change = config["fold_change"]
    log:
        TUXEDO_DIR + "03-flag/.log/{pheno}_tuxedo_flag.log"
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_flag.py "
        " -f {params.fold_change} "
        " {input.cuffdiff_gene_exp} "
        " {output.gene_flagged} "
        " 2>&1 | tee {log} && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_flag.py "
        " -f {params.fold_change} "
        " {input.cuffdiff_isoform_exp} "
        " {output.isoform_flagged} "
        " 2>&1 | tee >>{log} "
