rule ballgown_plots:
    input:
        BALLGOWN_DIR + '01-ballgown_diffex/ballgown_data.rda',
    output:
        expand(BALLGOWN_DIR + '01-ballgown_diffex/plots/by_phenotype/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf', phenotype = PHENOTYPES, dim = ['12','23'], ngenes = ['100','500']),
        expand(BALLGOWN_DIR + '01-ballgown_diffex/plots/by_phenotype/{phenotype}/MDSplot_{dim}_top{ngenes}.pdf', phenotype = PHENOTYPES, dim = ['12','23'], ngenes = ['100','500']),
        expand(BALLGOWN_DIR + '01-ballgown_diffex/plots/by_phenotype/{phenotype}/ScreePlot_top{ngenes}.pdf', phenotype = PHENOTYPES, ngenes = ['100','500']),
        expand(BALLGOWN_DIR + '01-ballgown_diffex/plots/by_phenotype/{phenotype}/BoxPlot.pdf', phenotype = PHENOTYPES),
        expand(BALLGOWN_DIR + '01-ballgown_diffex/plots/by_phenotype/{phenotype}/SampleHeatmap.pdf', phenotype = PHENOTYPES),
        expand(BALLGOWN_DIR + '01-ballgown_diffex/plots/by_phenotype/{phenotype}/Heatmap_TopVar.pdf', phenotype = PHENOTYPES),
        expand(BALLGOWN_DIR + '01-ballgown_diffex/plots/by_phenotype/{phenotype}/Heatmap_TopExp.pdf', phenotype = PHENOTYPES),
        expand(BALLGOWN_DIR + '01-ballgown_diffex/plots/comparison_plots/{phenotype}/VolcanoPlot_{comparison}.pdf',
               zip,
               phenotype=REPLICATE_PHENOTYPE_NAMES,
               comparison=REPLICATE_COMPARISON_GROUPS),
    log:
        BALLGOWN_DIR + '01-ballgown_diffex/.log/ballgown_plots.log'
    conda:
        'envs/ballgown_diffex.yaml'
    params:
        configfile_path = CONFIGFILE_PATH
    shell:
        '''
        Rscript {WATERMELON_SCRIPTS_DIR}/diffex_plots.R \
            --diffex_rda {input} \
            --config_file {params.configfile_path}
        rm -f Rplots.pdf #Some part of R generates this empty (nuisance) plot
        '''
