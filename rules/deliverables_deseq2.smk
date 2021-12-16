rule deliverables_deseq2:
    input:
        #counts
        counts = expand(DIFFEX_DIR + 'counts/deseq2_{name}.txt',
            name=['raw_counts', 'depth_normalized_counts', 'rlog_normalized_counts']),
        #Gene lists
        annot = helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + 'diffex_{model_name}/{contrast}.annot.txt',
            DESEQ2_CONTRAST_DICT),
        excel = helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + "diffex_{model_name}/{contrast}.annot.xlsx",
            DESEQ2_CONTRAST_DICT),
        #plots
        pca = expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf',
                phenotype = PHENOTYPES,
                dim = ['12','23'],
                ngenes = ['100','500']),
        scree = expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/ScreePlot_top{ngenes}.pdf',
                phenotype = PHENOTYPES,
                ngenes = ['100','500']),
        boxplots = expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/{plot_type}.pdf',
            phenotype = PHENOTYPES,
            plot_type = ['BoxPlot_raw', 'BoxPlot_rlog']),
        summaryplots = expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{plot_type}.pdf',
            plot_type = ['SampleHeatmap', 'Heatmap_TopVar', 'Heatmap_TopExp']),
        volcanoplots = helper.expand_model_contrast_filenames(\
                DIFFEX_DIR + 'diffex_{model_name}/volcano_plots/VolcanoPlot_{contrast}.pdf',
                DESEQ2_CONTRAST_DICT),
        summary_txt = DIFFEX_DIR + "summary/deseq2_summary.txt",
        summary_xlsx = DIFFEX_DIR + "summary/deseq2_summary.xlsx",
    output:
        #counts
        counts = expand(DELIVERABLES_DIR + 'counts/deseq2_{name}.txt',
            name=['raw_counts', 'depth_normalized_counts', 'rlog_normalized_counts']),
        #Gene lists - annotated and excel are placed in the same location in deliverables output
        annot = helper.expand_model_contrast_filenames(\
            DELIVERABLES_DIR + 'diffex_{model_name}/{contrast}.annot.txt',
            DESEQ2_CONTRAST_DICT),
        excel = helper.expand_model_contrast_filenames(\
            DELIVERABLES_DIR + "diffex_{model_name}/{contrast}.annot.xlsx",
            DESEQ2_CONTRAST_DICT),
        #plots
        pca = expand(DELIVERABLES_DIR + 'plots_labeled_by_pheno/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf',
                phenotype = PHENOTYPES,
                dim = ['12','23'],
                ngenes = ['100','500']),
        scree = expand(DELIVERABLES_DIR + 'plots_labeled_by_pheno/{phenotype}/ScreePlot_top{ngenes}.pdf',
                phenotype = PHENOTYPES,
                ngenes = ['100','500']),
        boxplots = expand(DELIVERABLES_DIR + 'plots_labeled_by_pheno/{phenotype}/{plot_type}.pdf',
            phenotype = PHENOTYPES,
            plot_type = ['BoxPlot_raw', 'BoxPlot_rlog']),
        summaryplots = expand(DELIVERABLES_DIR + 'plots_labeled_by_pheno/{plot_type}.pdf',
            plot_type = ['SampleHeatmap', 'Heatmap_TopVar', 'Heatmap_TopExp']),
        volcanoplots = helper.expand_model_contrast_filenames(\
                DELIVERABLES_DIR + 'diffex_{model_name}/volcano_plots/VolcanoPlot_{contrast}.pdf',
                DESEQ2_CONTRAST_DICT),
        summary_txt = DELIVERABLES_DIR + "summary/deseq2_summary.txt",
        summary_xlsx = DELIVERABLES_DIR + "summary/deseq2_summary.xlsx",
#    log:
#        DELIVERABLES_DIR + ".deliverables_deseq2.log"
    params:
        project_name = config['report_info']['project_name'],
        deseq2_input_dir = DIFFEX_DIR,
        deseq2_output_dir = DELIVERABLES_DIR,
        #Counts are placed in a separate location from the other deseq2 results
        counts_output_dir = DELIVERABLES_DIR + "counts",
    shell:
        """rsync -rlpgoD --include "diffex_*/*.annot.txt" --exclude counts --exclude "diffex_*/*.txt" --exclude ".*" --exclude "*.rda" {params.deseq2_input_dir} {params.deseq2_output_dir}
        for i in {input.counts} ; do cp $i {params.counts_output_dir} ; done"""
