rule tuxedo_cummerbund:
    input:
        genome_checksum = CONFIG_CHECKSUMS_DIR + "config-genome.watermelon.md5",
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        group_replicates = TUXEDO_DIR + "05-group_replicates/{pheno}/group_replicates.txt",
        gtf_file = "references/gtf"
    output:
        TUXEDO_DIR + "06-cummerbund/{pheno}/Plots",
        TUXEDO_DIR + "06-cummerbund/{pheno}/Plots/{pheno}_boxplot.pdf",
        TUXEDO_DIR + "06-cummerbund/{pheno}/{pheno}_repRawCounts.txt"
    params:
        cuff_diff_dir = TUXEDO_DIR + "01-cuffdiff/{pheno}",
        output_dir = TUXEDO_DIR + "06-cummerbund/{pheno}",
        genome = config["genome"],
    log:
         TUXEDO_DIR + "06-cummerbund/.log/{pheno}_cummerbund.log"
    shell:
        "module purge && module load watermelon_dependencies && "
        "mkdir -p {params.output_dir}/Plots && "
        "Rscript {WATERMELON_SCRIPTS_DIR}/Run_cummeRbund.R "
        " baseDir={params.output_dir} "
        " cuffDiffDir={params.cuff_diff_dir} "
        " grpRepFile={input.group_replicates} "
        " gtfFile={input.gtf_file} "
        " genome={params.genome} "
        " 2>&1 | tee {log} && "
        "touch {params.output_dir}/Plots "
