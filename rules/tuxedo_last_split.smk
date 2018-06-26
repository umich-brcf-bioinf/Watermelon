rule tuxedo_last_split:
    input:
        expand(TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_gene.txt", \
               zip, \
               phenotype_name=ALL_PHENOTYPE_NAMES, \
               comparison=ALL_COMPARISON_GROUPS)
    output:
        touch(TUXEDO_DIR + "07-split/last_split")
