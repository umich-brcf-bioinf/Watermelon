library(optparse)

option_list = list(
    make_option('--diffex_rda', type='character', help='[Required] Path to diffex rdata'),
    make_option('--config_file', type='character', help='[Required] Path to configfile')
)
opt = parse_args(OptionParser(option_list=option_list))

diffex_rda = opt$diffex_rda
config_file = opt$config_file

########################################################

library(ballgown)
library(dplyr)
library(GGally)
library(ggfortify)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(yaml)

########################################################
# Define plotting functions


plot_boxplot = function(mat, title, y_label, out_name = 'BoxPlot.pdf') {
    tidy_mat = tidyr::gather(as_tibble(mat), key = 'sample', value = 'counts')

    box_plot = ggplot(tidy_mat, aes(x = sample, y = counts)) +
        geom_boxplot(notch = TRUE) +
        labs(
            title = title,
            x = '',
            y = y_label) +
        theme_bw() + theme(axis.text.x = element_text(angle = 90))
    ggsave(filename = file.path(plots_dir, 'summary_plots', out_name), plot = box_plot, height = 8, width = 8, dpi = 300)

    return(box_plot)
}


plot_sample_correlation_heatmap = function(mat, out_name = 'SampleHeatmap.pdf') {

    dist_obj = dist(t(mat), method = 'euclidean', diag = T, upper = T)
    hclust_obj = hclust(dist_obj, method = 'complete')
    dist_mat = as.matrix(dist_obj)

    pdf(file = file.path(plots_dir, 'summary_plots', out_name), height = 8, width = 8)
        pheatmap(
            mat = dist_mat,
            cluster_rows = hclust_obj,
            cluster_cols = hclust_obj,
            color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Blues')))(255)
        )
    dev.off()

    return(list(dist_obj = dist_obj, hclust_obj = hclust_obj))
}


plot_top_variably_expressed_heatmap = function() {

}


plot_top_expressed_heatmap = function() {

}


add_replicate_col = function(pdata, factor_name) {
    if(!(factor_name %in% colnames(pdata))) {
        stop(sprintf('factor_name = %s is not in columns of pdata.', factor_name))
    }

    # Add replicate column for the factor_name
    pdata_list = lapply(split(pdata, pdata[, factor_name]), function(l) {
        l$replicate = seq_along(l$sample)
        return(l)
    })
    pdata = Reduce(rbind, pdata_list)
    pdata$replicate = factor(pdata$replicate)

    return(pdata)
}


plot_PCA = function(mat, pdata, factor_name, top_n = 1000, out_name = 'PCAplot.pdf') {

    # Add replicate column
    pdata = add_replicate_col(pdata, factor_name)

    # Calculate PCA on top_n variable genes
    top_var_mat = mat[order(matrixStats::rowVars(mat), decreasing = T), ][1:top_n, ]

    pca = prcomp(top_var_mat, retx = TRUE, scale = TRUE)
    pca_df = data.frame(
        sample = rownames(pca$rotation),
        PC1 = pca$rotation[, 'PC1'],
        PC2 = pca$rotation[, 'PC2'],
        stringsAsFactors = F
    )

    # Compute variance explained
    var_explained = round((pca$sdev^2 / sum(pca$sdev^2))*100, 1)

    # Join the pca_df with pdata to get the replicates and factor_name columns
    pca_df = merge(pca_df, pdata, by = 'sample', all.x = T)

    # Plot PCA
    pca_plot = ggplot(pca_df, aes(x = PC1, y = PC2, color = replicate)) +
        geom_point(aes_string(shape = factor_name)) +
        labs(
            title = sprintf('%s PCA plot', factor_name),
            subtitle = sprintf('Using top %s variable genes', top_n),
            x = sprintf('PC1: %s%% variance', var_explained[1]),
            y = sprintf('PC2: %s%% variance', var_explained[2])) +
        theme_bw()
    ggsave(filename = file.path(plots_dir, 'comparison_plots', factor_name, out_name), plot = pca_plot, height = 8, width = 8, dpi = 300)

    return(pca_plot)
}


plot_MDS = function(mat, pdata, factor_name, top_n = 1000, out_name = 'MDSplot.pdf') {

    # Add replicate column
    pdata = add_replicate_col(pdata, factor_name)

    # Calculate MDS and pull relevant parts for the plot
    mds = limma::plotMDS(x = mat, top = top_n, plot = FALSE, gene.selection = 'common')
    mds_df = data.frame(
        sample = rownames(mds$cmdscale.out),
        x = mds$x,
        y = mds$y,
        stringsAsFactors = F
    )

    # Join the mds_df with pdata to get the replicates and factor_name columns
    mds_df = merge(mds_df, pdata, by = 'sample', all.x = T)

    # Plot MDS
    mds_plot = ggplot(mds_df, aes(x = x, y = y, color = replicate)) +
        geom_point(aes_string(shape = factor_name)) +
        labs(
            title = sprintf('%s MDS plot', factor_name),
            subtitle = sprintf('Using top %s variable genes', top_n),
            x = 'Dim 1',
            y = 'Dim 2') +
        theme_bw()
    ggsave(filename = file.path(plots_dir, 'comparison_plots', factor_name, out_name), plot = mds_plot, height = 8, width = 8, dpi = 300)

    return(mds_plot)
}


plot_volcano = function(de_list, method = c('ballgown', 'deseq2'), comparison_type, exp_name, con_name, fdr_cutoff, logfc_cutoff, out_name) {

    method = match.arg(method)

    # Determine the correct columns on the basis of deseq2 or ballgown
    if(method == 'ballgown') {
        log2fc = 'log2fc'
        pval = 'pval'
        padj = 'qval'
        gene_id = 'id'
    } else {
        stop('Sorry, deseq2 is not yet supported.')
    }

    # Add direction column
    de_list$direction = 'NS'
    de_list$direction[de_list$diff_exp == 'YES' & de_list$log2fc <= 0] = 'Down'
    de_list$direction[de_list$diff_exp == 'YES' & de_list$log2fc > 0] = 'Up'
    de_list$direction = factor(de_list$direction, levels = c('Up', 'Down', 'NS'))

    if(unique(de_list$direction == 'NS')) {
        warning(sprintf('No genes were DE at fdr < %s and |log2fc| > %s. Consider a different threshold.', fdr_cutoff, logfc_cutoff))
    }

    # Transform qval to -log10 scale
    de_list$log10qval = -log10(de_list$qval)

    # Add top 10 Up and 10 Down gene labels
    top = rbind(
        head(subset(de_list, direction == 'Up'), 10),
        head(subset(de_list, direction == 'Down'), 10))
    top$label = top$id
    de_list = merge(x = de_list, y = top[, c('id','label')], by = 'id', all.x = TRUE, sort = FALSE)

    # Volcano Plot
    volcano_plot = ggplot(de_list, aes_string(x = 'log2fc', y = 'log10qval', color = 'direction')) +
        geom_point(size = 1) +
        scale_color_manual(name = '', values=c('#B31B21', '#1465AC', 'darkgray')) +
        geom_vline(xintercept = c(0, -1*fc_cutoff, fc_cutoff), linetype = c(1, 2, 2), color = c('black', 'black', 'black')) +
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

    ggsave(filename = file.path(plots_dir, 'comparison_plots', comparison_type, out_name), plot = volcano_plot)

    return(volcano_plot)
}



########################################################
########################################################
########################################################
# Parse yaml file

yaml = read_yaml(config_file)

fdr_cutoff = yaml$deseq2_adjustedPValue
logfc_cutoff = log2(yaml$fold_change)

diffex_dir = yaml$diffex_output_dir
results_dir = sprintf('%s/ballgown/01-ballgown_diffex', diffex_dir)
plots_dir = sprintf('%s/plots', results_dir)

########################################################
# Load data

load(diffex_rda)

if('bg_data' %in% ls()) {
    method = 'ballgown'
    boxplot_y_lab = 'log2(FPKM)'
    pdata = pData(bg_data)
} else {
    method = 'deseq2'
    boxplot_y_lab = 'log2(counts)'
}

# Setup FPKM matrix for use in PCA and heatmaps
if(method == 'ballgown') {
    mat = log2(as.matrix(gene_fpkms[,-1]) + 1)
    colnames(mat) = gsub('FPKM.','', colnames(mat))
} else {

}

########################################################
# Plots

# Boxplot

log2_boxplot = plot_boxplot(mat = mat, title = 'log2 FPKMs', y_label = 'log2(FPKM)')

########################################################
# Heatmaps

log2_heatmap = plot_sample_correlation_heatmap(mat = mat)

########################################################
# PCA and MDS plots for each of the factor_names

factor_names = setdiff(colnames(pData(bg_data)), 'sample')
for(factor_name in factor_names) {

    log2_pca = plot_PCA(mat = mat, pdata = pdata, factor_name = factor_name, top_n = 1000, out_name = 'PCAplot.pdf')
    log2_mds = plot_MDS(mat = mat, pdata = pdata, factor_name = factor_name, top_n = 1000, out_name = 'MDSplot.pdf')

}

########################################################
# Volcano plot for each comparison

for(comparison_type in comparison_types) {
    # What are all the comparisons to be done based on the type?
    all_comparisons = comparisons[[comparison_type]]

    # Do the gene-level and isoform-level tests for DE
    for(i in 1:nrow(all_comparisons)) {
        # Extract experiment name and control name for convenience
        exp = as.character(all_comparisons[i, 'exp'])
        con = as.character(all_comparisons[i, 'con'])
        out_name = as.character(all_comparisons[i, 'out_name'])

        message(sprintf('Plotting volcano plot for %s %s', comparison_type, out_name))

        gene_list = de_results[[comparison_type]][[out_name]][['gene_list']]

        volcano_plot = plot_volcano(
            de_list = gene_list,
            method = method,
            comparison_type = comparison_type,
            exp_name = exp,
            con_name = con,
            fdr_cutoff = fdr_cutoff,
            logfc_cutoff = logfc_cutoff,
            out_name = sprintf('VolcanoPlot_%s.pdf', out_name))
    }
}
