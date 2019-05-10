rule ballgown_excel:
    input:
        gene_annot = BALLGOWN_DIR + '02-annotate/{phenotype}/{comparison}_gene.annot.txt',
        isoform_annot = BALLGOWN_DIR + '02-annotate/{phenotype}/{comparison}_isoform.annot.txt',
        glossary = BALLGOWN_DIR + '03-run_info/glossary.txt',
        run_info = BALLGOWN_DIR + '03-run_info/run_info.txt'
    output:
        BALLGOWN_DIR + '04-excel/{phenotype}/{comparison}.xlsx'
    log:
        BALLGOWN_DIR + '04-excel/.log/{phenotype}_{comparison}.ballgown_excel.log'
    shell:
        '''
        module purge && module load python/3.6.1
        python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py \
         -g {input.gene_annot} \
         -i {input.isoform_annot} \
         --glossary {input.glossary} \
         --info_filepath {input.run_info} \
         {output} \
         2>&1 | tee {log}
        '''
