rule ballgown_plots:
    input:
        BALLGOWN_DIR + '01-ballgown_diffex/ballgown_data.rda',
    output:
        expand(BALLGOWN_DIR + '01-ballgown_diffex/plots/comparison_plots/{phenotype}/PCAplot.pdf', phenotype = PHENOTYPES),
        expand(BALLGOWN_DIR + '01-ballgown_diffex/plots/comparison_plots/{phenotype}/MDSplot.pdf', phenotype = PHENOTYPES),
        BALLGOWN_DIR + '01-ballgown_diffex/plots/summary_plots/BoxPlot.pdf',
        BALLGOWN_DIR + '01-ballgown_diffex/plots/summary_plots/SampleHeatmap.pdf',
        BALLGOWN_DIR + '01-ballgown_diffex/plots/summary_plots/Heatmap_TopVar.pdf',
        BALLGOWN_DIR + '01-ballgown_diffex/plots/summary_plots/Heatmap_TopExp.pdf',
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
