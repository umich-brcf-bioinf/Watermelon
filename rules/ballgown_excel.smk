rule ballgown_excel:
    input:
        gene_annot = BALLGOWN_DIR + '02-annotate/{phenotype}/{comparison}_gene.annot.txt',
        isoform_annot = BALLGOWN_DIR + '02-annotate/{phenotype}/{comparison}_isoform.annot.txt',
        glossary = WATERMELON_SCRIPTS_DIR + 'ballgown_glossary.txt',
    output:
        BALLGOWN_DIR + '04-excel/{phenotype}/{comparison}.xlsx'
    log:
        BALLGOWN_DIR + '04-excel/.log/{phenotype}_{comparison}.ballgown_excel.log'
    conda:
        '../envs/python_3.6.1.yaml'
    shell:
        '''
        python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py \
         -g {input.gene_annot} \
         -i {input.isoform_annot} \
         --glossary {input.glossary} \
         {output} \
         2>&1 | tee {log}
        '''
