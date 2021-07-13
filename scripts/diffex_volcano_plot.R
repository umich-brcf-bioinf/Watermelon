##########
# Set up logging and save snakemake S4 object (for debugging or running manually)
log = file(snakemake@log[[1]], open = 'wt')
sink(log)
sink(log, type = 'message')
save(snakemake, file = snakemake@params[['snakemake_rdata']])

#setwd("/nfs/med-bfx-activeprojects/trsaari/sandbox/20190718_test_lumeng/")
#load("/nfs/med-bfx-activeprojects/trsaari/sandbox/20190718_test_lumeng/analysis_20190920/diffex_results/deseq2/plots/comparison_plots/model_one/.comparison_plot_DM.fem_v_NDM.fem_snakemake.rda")

#Isolate conda environment: https://github.com/conda-forge/r-base-feedstock/issues/37
#If we move away from conda in the future, we may want to remove this
.libPaths(R.home("library"))

########################################################
# Load libraries
lib.vector = c("data.table", "DESeq2", "dplyr", "GGally", "ggfortify", "ggplot2",
    "ggrepel", "pheatmap", "RColorBrewer", "stringr", "tidyr", "yaml")

foo = suppressMessages(lapply(lib.vector, library, character.only=T, warn.conflicts=F, quietly=T))

########################################################
# Functions
plot_volcano = function(de_list, method = c('ballgown', 'deseq2'), exp_name, con_name, fdr_cutoff, logfc_cutoff, out_filepath_pdf, out_filepath_png) {

    method = match.arg(method)

    # Determine the correct column names on the basis of deseq2 or ballgown
    if(method == 'deseq2') {
        log2fc = 'log2FoldChange'
        pval = 'pvalue'
        padj = 'padj'
        id = 'gene_id'
        de_call = 'Call'
    } else {
        log2fc = 'log2fc'
        pval = 'pval'
        padj = 'qval'
        id = 'gene_id'
        de_call = 'diff_exp'
    }

    # Add direction column
    de_list$direction = 'NS'
    de_list$direction[de_list[, de_call] == 'YES' & de_list[, log2fc] <= 0] = 'Down'
    de_list$direction[de_list[, de_call] == 'YES' & de_list[, log2fc] > 0] = 'Up'

    # Determine direction labels (with number per)
    direction_table = table(de_list$direction)

    if('Up' %in% names(direction_table)) {
        up_label = sprintf('Up: %s', direction_table['Up'])
    } else {
        up_label = 'Up: 0'
    }

    if('Down' %in% names(direction_table)) {
        down_label = sprintf('Down: %s', direction_table['Down'])
    } else {
        down_label = 'Down: 0'
    }

    if('NS' %in% names(direction_table)) {
        ns_label = sprintf('NS: %s', direction_table['NS'])
    } else {
        ns_label = 'NS: 0'
    }

    de_list$direction_count = factor(
        de_list$direction,
        levels = c('Up', 'Down', 'NS'),
        labels = c(up_label, down_label, ns_label))

    if(all( !( c('Up', 'Downl') %in% names(direction_table) ) )) {
        warning(sprintf('No genes were DE at fdr < %s and |log2fc| > %s. Consider a different threshold.', fdr_cutoff, logfc_cutoff))
    }

    # Transform qval to -log10 scale
    de_list$log10qval = -log10(de_list[, padj])

    # Add top 10 Up and 10 Down gene labels
    # de_list is assumed to be ordered by Call/diff_exp and then qvalue from deseq2_diffex.R and ballgown_diffex.R
    top = rbind(
        head(subset(de_list, direction == 'Up'), 10),
        head(subset(de_list, direction == 'Down'), 10))
    top$label = top$external_gene_name
    de_list = merge(x = de_list, y = top[, c(id,'label')], by = id, all.x = TRUE, sort = FALSE)

    # Volcano Plot
    volcano_plot = ggplot(de_list, aes_string(x = log2fc, y = 'log10qval', color = 'direction_count', alpha = 0.5)) +
        scale_alpha(guide = 'none') +
        geom_point(size = 1) +
        scale_color_manual(name = '', values=c('#B31B21', '#1465AC', 'darkgray')) +
        geom_vline(xintercept = c(0, -1*logfc_cutoff, logfc_cutoff), linetype = c(1, 2, 2), color = c('black', 'black', 'black')) +
        geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, color = 'black') +
        labs(
            title = sprintf('%s_v_%s', exp_name, con_name),
            x = 'Log2 fold-change',
            y = '-Log10 adjusted q-value') +
        theme_classic()
    # Add gene symbol labels
    if(!all(is.na(de_list$label))) {
        volcano_plot = volcano_plot + ggrepel::geom_label_repel(label = de_list$label, force = 3, segment.alpha = 0.4)
    }
    ggsave(filename = out_filepath_pdf, plot = volcano_plot, height = 8, width = 8, dpi = 300)
    ggsave(filename = out_filepath_png, plot = volcano_plot, height = 8, width = 8, dpi = 300)

    return(volcano_plot)
}

###########################
# Set up for volcano plots

# Load gene list
gene_list = "foo" # Needs to be the annotated table
#gene_list = data.frame(de_results[[comparison_type]][[comparison]][['gene_list']])

# Volcano plot
out_pdf = snakemake@output[['volcano_plot_pdf']]
out_png = snakemake@output[['volcano_plot_png']]


######################
# Create volcano plot
message(sprintf('Plotting volcano plot for %s %s', factor_name, contrast_name))

volcano_plot = plot_volcano(
  de_list = gene_list,
  method = "deseq2",
  exp_name = test_name,
  con_name = reference_name,
  fdr_cutoff = fdr_cutoff,
  logfc_cutoff = fc_cutoff,
  out_filepath_pdf = out_pdf,
  out_filepath_png = out_png)
