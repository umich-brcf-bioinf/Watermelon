rule deseq2_metadata_contrasts:
    input:
        sample_checksum      = CONFIG_CHECKSUMS_DIR + "config-samples.watermelon.md5",
        comparison_checksum  = CONFIG_CHECKSUMS_DIR + "config-comparisons.watermelon.md5",
        phenotype_checksum   = CONFIG_CHECKSUMS_DIR + "config-phenotypes.watermelon.md5",
        main_factor_checksum = CONFIG_CHECKSUMS_DIR + "config-main_factors.watermelon.md5"
    output:
        sample_metadata = DESEQ2_DIR + "01-metadata_contrasts/sample_metadata.txt",
        contrasts = DESEQ2_DIR + "01-metadata_contrasts/contrasts.txt"
    params:
        phenos_with_replicates = phenotypeManager.phenotypes_with_replicates
    run:
        deseq2_helper.build_sample_metadata(config, params.phenos_with_replicates, output.sample_metadata)
        deseq2_helper.build_contrasts(config, params.phenos_with_replicates, output.contrasts)
