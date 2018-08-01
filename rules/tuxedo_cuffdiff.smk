rule tuxedo_cuffdiff:
    input:
        sample_checksum = CONFIG_CHECKSUMS_DIR + "phenotype_samples-{pheno}.watermelon.md5",
        comparison_checksum = CONFIG_CHECKSUMS_DIR + "phenotype_comparisons-{pheno}.watermelon.md5",
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        fasta_file = "references/bowtie2_index/genome.fa",
        gtf_file = "references/gtf",
        bam_files = expand(ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_accepted_hits.bam",
                           sample=config[SAMPLES_KEY])
    output:
        TUXEDO_DIR + "01-cuffdiff/{pheno}/gene_exp.diff",
        TUXEDO_DIR + "01-cuffdiff/{pheno}/isoform_exp.diff",
        TUXEDO_DIR + "01-cuffdiff/{pheno}/read_groups.info"
    params:
        output_dir = TUXEDO_DIR + "01-cuffdiff/{pheno}",
        labels = lambda wildcards : phenotypeManager.concatenated_comparison_values(',')[wildcards.pheno],
        samples = lambda wildcards : phenotypeManager.cuffdiff_samples(wildcards.pheno,
                                                                       ALIGNMENT_DIR + "04-tophat/{sample_placeholder}/{sample_placeholder}_accepted_hits.bam"),
        strand = rnaseq_snakefile_helper.strand_option_tophat(config["alignment_options"]["library_type"])
    threads: 8
    log:
        TUXEDO_DIR + "01-cuffdiff/.log/{pheno}_cuffdiff.log"
    shell:
        "rm -rf {params.output_dir} {params.output_dir}.tmp &&"
        "module purge && module load watermelon_dependencies && "
        "cuffdiff -q "
        " -p {threads} "
        " -L {params.labels} "
        " --max-bundle-frags 999999999 "
        " --library-type {params.strand} "
        " -o {params.output_dir}.tmp "
        " -b {input.fasta_file} "
        " -u -N "
        " --compatible-hits-norm "
        " {input.gtf_file} "
        " {params.samples} "
        " 2>&1 | tee {log} && "
        "mv {params.output_dir}.tmp {params.output_dir} "
