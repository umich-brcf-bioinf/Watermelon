rule deliverables_deseq2:
    input:
        #counts
        counts = expand(DIFFEX_DIR + 'deseq2/counts/deseq2_{name}.txt',
            name=['raw_counts', 'depth_normalized_counts', 'rlog_normalized_counts']),
        #Gene lists
        gene_lists = rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + 'deseq2/gene_lists/{model_name}/{contrast}.txt',
            DESEQ2_CONTRAST_DICT),
        excel = rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + "deseq2/excel/{model_name}/{contrast}.xlsx",
            DESEQ2_CONTRAST_DICT),
        #plots
        pca = expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf',
                phenotype = PHENOTYPES,
                dim = ['12','23'],
                ngenes = ['100','500']),
        scree = expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/ScreePlot_top{ngenes}.pdf',
                phenotype = PHENOTYPES,
                ngenes = ['100','500']),
        summaryplots = expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/{plotType}.pdf',
            phenotype = PHENOTYPES, plotType = ['BoxPlot', 'SampleHeatmap', 'Heatmap_TopVar', 'Heatmap_TopExp']),
        volcanoplots = rnaseq_snakefile_helper.expand_model_contrast_filenames(\
                DIFFEX_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.pdf',
                DESEQ2_CONTRAST_DICT),
        summary_txt = DIFFEX_DIR + "deseq2/summary/deseq2_summary.txt",
        summary_xlsx = DIFFEX_DIR + "deseq2/summary/deseq2_summary.xlsx",
    output:
        #counts
        counts = expand(DELIVERABLES_DIR + 'counts/deseq2_{name}.txt',
            name=['raw_counts', 'depth_normalized_counts', 'rlog_normalized_counts']),
        #Gene lists
        gene_lists = rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DELIVERABLES_DIR + 'deseq2/gene_lists/{model_name}/{contrast}.txt',
            DESEQ2_CONTRAST_DICT),
        excel = rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DELIVERABLES_DIR + "deseq2/excel/{model_name}/{contrast}.xlsx",
            DESEQ2_CONTRAST_DICT),
        #plots
        pca = expand(DELIVERABLES_DIR + 'deseq2/plots/by_phenotype/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf',
                phenotype = PHENOTYPES,
                dim = ['12','23'],
                ngenes = ['100','500']),
        scree = expand(DELIVERABLES_DIR + 'deseq2/plots/by_phenotype/{phenotype}/ScreePlot_top{ngenes}.pdf',
                phenotype = PHENOTYPES,
                ngenes = ['100','500']),
        summaryplots = expand(DELIVERABLES_DIR + 'deseq2/plots/by_phenotype/{phenotype}/{plotType}.pdf',
            phenotype = PHENOTYPES, plotType = ['BoxPlot', 'SampleHeatmap', 'Heatmap_TopVar', 'Heatmap_TopExp']),
        volcanoplots = rnaseq_snakefile_helper.expand_model_contrast_filenames(\
                DELIVERABLES_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.pdf',
                DESEQ2_CONTRAST_DICT),
        summary_txt = DELIVERABLES_DIR + "deseq2/summary/deseq2_summary.txt",
        summary_xlsx = DELIVERABLES_DIR + "deseq2/summary/deseq2_summary.xlsx",
#    log:
#        DELIVERABLES_DIR + ".deliverables_deseq2.log"
    params:
        #Counts are placed in a separate location from the other deseq2 results
        deseq2_input_dir = DIFFEX_DIR + "deseq2/",
        deseq2_output_dir = DELIVERABLES_DIR + "deseq2",
        counts_output_dir = DELIVERABLES_DIR + "counts"
    shell:
        """rsync -rlpgoD --exclude counts --exclude ".*" --exclude "*.rda" --exclude "*.annot.txt" {params.deseq2_input_dir} {params.deseq2_output_dir}
        for i in {input.counts} ; do cp $i {params.counts_output_dir} ; done"""
