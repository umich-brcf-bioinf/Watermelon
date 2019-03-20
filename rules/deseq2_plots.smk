rule deseq2_plots:
    input:
        DESEQ2_DIR + '02-deseq2_diffex/deseq2_data.rda',
    output:
        expand(DESEQ2_DIR + '02-deseq2_diffex/plots/by_phenotype/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf', phenotype = PHENOTYPES, dim = ['12','23'], ngenes = ['100','500']),
        expand(DESEQ2_DIR + '02-deseq2_diffex/plots/by_phenotype/{phenotype}/MDSplot_{dim}_top{ngenes}.pdf', phenotype = PHENOTYPES, dim = ['12','23'], ngenes = ['100','500']),
        expand(DESEQ2_DIR + '02-deseq2_diffex/plots/by_phenotype/{phenotype}/ScreePlot_top{ngenes}.pdf', phenotype = PHENOTYPES, ngenes = ['100','500']),
        expand(DESEQ2_DIR + '02-deseq2_diffex/plots/by_phenotype/{phenotype}/BoxPlot.pdf', phenotype = PHENOTYPES),
        expand(DESEQ2_DIR + '02-deseq2_diffex/plots/by_phenotype/{phenotype}/SampleHeatmap.pdf', phenotype = PHENOTYPES),
        expand(DESEQ2_DIR + '02-deseq2_diffex/plots/by_phenotype/{phenotype}/Heatmap_TopVar.pdf', phenotype = PHENOTYPES),
        expand(DESEQ2_DIR + '02-deseq2_diffex/plots/by_phenotype/{phenotype}/Heatmap_TopExp.pdf', phenotype = PHENOTYPES),
        expand(DESEQ2_DIR + '02-deseq2_diffex/plots/comparison_plots/{phenotype}/VolcanoPlot_{comparison}.pdf',
               zip,
               phenotype=REPLICATE_PHENOTYPE_NAMES,
               comparison=REPLICATE_COMPARISON_GROUPS),
    log:
        DESEQ2_DIR + '02-deseq2_diffex/.log/deseq2_plots.log'
    conda:
        'envs/diffex.yaml'
    params:
        configfile_path = CONFIGFILE_PATH
    shell:
        '''(Rscript {WATERMELON_SCRIPTS_DIR}/diffex_plots.R \
            --diffex_rda {input} \
            --config_file {params.configfile_path}
        rm -f Rplots.pdf #Some part of R generates this empty (nuisance) plot
        ) 2>&1 | tee {log}'''
