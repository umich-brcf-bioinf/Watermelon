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
library(DESeq2)
library(dplyr)
library(GGally)
library(ggfortify)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(tidyr)
library(yaml)

########################################################
# Define plotting functions


plot_boxplot = function(mat, pdata, factor_name, title, y_label, out_name = 'BoxPlot.pdf') {

    annot_df = data.frame(
        sample = pdata$sample,
        Group = factor(pdata[, factor_name]),
        row.names = pdata$sample,
        stringsAsFactors = F
    )

    tidy_mat = tidyr::gather(as_tibble(mat), key = 'sample', value = 'counts') %>%
        left_join(annot_df, by = 'sample')

    box_plot = ggplot(tidy_mat, aes(x = sample, y = counts, fill = Group)) +
        geom_boxplot(notch = TRUE) +
        labs(
            title = title,
            x = '',
            y = y_label) +
        theme_bw() + theme(axis.text.x = element_text(angle = 90))
    ggsave(filename = file.path(plots_dir, 'by_phenotype', factor_name, out_name), plot = box_plot, height = 8, width = 8, dpi = 300)

    return(box_plot)
}


plot_sample_correlation_heatmap = function(mat, pdata, factor_name, out_name = 'SampleHeatmap.pdf') {

    if(!(factor_name %in% colnames(pdata))) {
        stop(sprintf('factor_name = %s is not in columns of pdata.', factor_name))
    }

    annot_df = data.frame(
        Group = pdata[, factor_name],
        row.names = pdata$sample
    )

    # Calculate the distance matrix and clustering object
    dist_obj = dist(t(mat), method = 'euclidean', diag = T, upper = T)
    hclust_obj = hclust(dist_obj, method = 'complete')
    dist_mat = as.matrix(dist_obj)

    pdf(file = file.path(plots_dir, 'by_phenotype', factor_name, out_name), height = 8, width = 8)
        pheatmap(
            mat = dist_mat,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            annotation_col = annot_df,
            color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Blues')))(255)
        )
    dev.off()

    return(list(dist_obj = dist_obj, hclust_obj = hclust_obj))
}


plot_top_variably_expressed_heatmap = function(mat, pdata, factor_name, top_n = 1000, out_name = 'Heatmap_TopVar.pdf') {

    if(!(factor_name %in% colnames(pdata))) {
        stop(sprintf('factor_name = %s is not in columns of pdata.', factor_name))
    }

    annot_df = data.frame(
        Group = pdata[, factor_name],
        row.names = pdata$sample
    )

    # Calculate the top_n variable genes and put them in decreasing order
    top_var_mat = mat[order(matrixStats::rowVars(mat), decreasing = T), ][1:top_n, ]

    pdf(file = file.path(plots_dir, 'by_phenotype', factor_name, out_name), height = 20, width = 10)
        pheatmap(
            mat = top_var_mat,
            scale= 'row',
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = FALSE,
            annotation_col = annot_df,
            fontsize = 7,
            fontsize_row = 7,
            las = 2,
            main = sprintf('Top %s variably expressed genes', top_n),
            color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Blues')))(255)
        )
    dev.off()
}


plot_top_expressed_heatmap = function(mat, pdata, factor_name, top_n = 1000, out_name = 'Heatmap_TopExp.pdf') {

    if(!(factor_name %in% colnames(pdata))) {
        stop(sprintf('factor_name = %s is not in columns of pdata.', factor_name))
    }

    annot_df = data.frame(
        Group = pdata[, factor_name],
        row.names = pdata$sample
    )

    # Calculate the top_n expressed genes and put them in decreasing order
    top_exp_mat = mat[order(rowMeans(mat), decreasing = T), ][1:top_n, ]

    pdf(file = file.path(plots_dir, 'by_phenotype', factor_name, out_name), height = 20, width = 10)
        pheatmap(
            mat = top_exp_mat,
            scale= 'row',
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = FALSE,
            annotation_col = annot_df,
            fontsize = 7,
            fontsize_row = 7,
            las = 2,
            main = sprintf('Top %s expressed genes', top_n),
            color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Blues')))(255)
        )
    dev.off()
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


compute_PCA = function(mat, pdata, factor_name, top_n = 500, dims = c('PC1','PC2')) {
    # Check for validly formatted dims
    if(!all(grepl('PC\\d{1,}', dims))) {
        stop("dims parameter must be e.g. c('PC1','PC2'). That is, any valid PC from prcomp.")
    }

    # Pull digits from dims
    pc_digits = as.integer(gsub('PC','',dims))

    # Add replicate column
    pdata = add_replicate_col(pdata, factor_name)

    # Calculate the top_n variable genes and put them in decreasing order
    top_var_mat = mat[order(matrixStats::rowVars(mat), decreasing = T), ][1:top_n, ]

    # Compute the PCA, I looked at the DESeq2 code to match the result in the previous plots
    pca = prcomp(t(top_var_mat), retx = TRUE)
    pca_df = data.frame(
        sample = rownames(pca$x),
        x = pca$x[, dims[1]],
        y = pca$x[, dims[2]],
        stringsAsFactors = F
    )

    # Compute variance explained and extract for the correct dimensions
    all_var_explained = round((pca$sdev^2 / sum(pca$sdev^2))*100, 1)
    var_explained = round((pca$sdev^2 / sum(pca$sdev^2))*100, 1)[pc_digits]

    # Join the pca_df with pdata to get the replicates and factor_name columns
    pca_df = merge(pca_df, pdata, by = 'sample', all.x = T)

    return(list(all_var_explained = all_var_explained, var_explained = var_explained, pca_df = pca_df, factor_name = factor_name, top_n = top_n, dims = dims))
}


plot_scree = function(compute_PCA_result, out_name = 'ScreePlot.pdf') {
    pca_df = compute_PCA_result[['pca_df']]
    all_var_explained = compute_PCA_result[['all_var_explained']]
    top_n = compute_PCA_result[['top_n']]
    factor_name = compute_PCA_result[['factor_name']]

    scree_tbl = as_tibble(data.frame(
        PC = factor(paste0('PC', seq_along(all_var_explained)), level = paste0('PC', seq_along(all_var_explained))),
        var_explained = all_var_explained, stringsAsFactors = F))
    scree_tbl = scree_tbl[1:min(10, nrow(scree_tbl)), ]

    # Plot scree
    scree_plot = ggplot(scree_tbl, aes(x = PC, y = var_explained)) + geom_bar(stat = 'identity') +
        labs(
            title = sprintf('%s Scree plot', factor_name),
            subtitle = sprintf('PCA computed with top %s variable genes', top_n),
            x = 'Principal Component',
            y = 'Percent Variance Explained'
        ) +
        theme_bw()
    ggsave(filename = file.path(plots_dir, 'by_phenotype', factor_name, out_name), plot = scree_plot, height = 6, width = 6, dpi = 300)

    return(scree_plot)
}


plot_PCA = function(compute_PCA_result, out_name = 'PCAplot.pdf') {

    pca_df = compute_PCA_result[['pca_df']]
    var_explained = compute_PCA_result[['var_explained']]
    top_n = compute_PCA_result[['top_n']]
    factor_name = compute_PCA_result[['factor_name']]
    dims = compute_PCA_result[['dims']]

    # Plot PCA
    pca_plot = ggplot(pca_df, aes(x = x, y = y, color = replicate)) +
        geom_point(aes_string(shape = factor_name)) +
        labs(
            title = sprintf('%s PCA plot', factor_name),
            subtitle = sprintf('Using top %s variable genes', top_n),
            x = sprintf('%s: %s%% variance', dims[1], var_explained[1]),
            y = sprintf('%s: %s%% variance', dims[2], var_explained[2])) +
        theme_bw()
    ggsave(filename = file.path(plots_dir, 'by_phenotype', factor_name, out_name), plot = pca_plot, height = 6, width = 6, dpi = 300)

    return(pca_plot)
}


plot_MDS = function(mat, pdata, factor_name, top_n = 500, dims = c(1,2), out_name = 'MDSplot.pdf') {

    if(!is(dims, 'numeric')) {
        stop('dims must be numeric, e.g. c(1,2).')
    }
    if(!all(dims < 10)) {
        stop('dims should not exceed the 9th dimension.')
    }

    # Add replicate column
    pdata = add_replicate_col(pdata, factor_name)

    # Calculate MDS and pull relevant parts for the plot
    mds = limma::plotMDS(x = mat, top = top_n, dim.plot = dims, plot = FALSE, gene.selection = 'common')
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
            x = sprintf('Dim %s', dims[1]),
            y = sprintf('Dim %s', dims[2])) +
        theme_bw()
    ggsave(filename = file.path(plots_dir, 'by_phenotype', factor_name, out_name), plot = mds_plot, height = 6, width = 6, dpi = 300)

    return(mds_plot)
}


plot_volcano = function(de_list, method = c('ballgown', 'deseq2'), comparison_type, exp_name, con_name, fdr_cutoff, logfc_cutoff, out_name) {

    method = match.arg(method)

    # Determine the correct column names on the basis of deseq2 or ballgown
    if(method == 'ballgown') {
        log2fc = 'log2fc'
        pval = 'pval'
        padj = 'qval'
        gene_id = 'id'
        de_call = 'diff_exp'
    } else {
        log2fc = 'log2FoldChange'
        pval = 'pvalue'
        padj = 'padj'
        gene_id = 'id'
        de_call = 'Call'
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
    top$label = top$id
    de_list = merge(x = de_list, y = top[, c('id','label')], by = 'id', all.x = TRUE, sort = FALSE)

    # Volcano Plot
    volcano_plot = ggplot(de_list, aes_string(x = log2fc, y = 'log10qval', color = 'direction_count')) +
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
    ggsave(filename = file.path(plots_dir, 'comparison_plots', comparison_type, out_name), plot = volcano_plot, height = 8, width = 8, dpi = 300)

    return(volcano_plot)
}



########################################################
########################################################
########################################################
# Parse yaml file

yaml = read_yaml(config_file)

fdr_cutoff = yaml$deseq2_adjustedPValue
logfc_cutoff = log2(yaml$fold_change)

phenotypes = str_trim(unlist(str_split(yaml$phenotypes, '\\^')))

diffex_dir = yaml$diffex_output_dir

########################################################
# Load data

load(diffex_rda)

if('bg_data' %in% ls()) {
    method = 'ballgown'

    results_dir = sprintf('%s/ballgown/01-ballgown_diffex', diffex_dir)
    plots_dir = sprintf('%s/plots', results_dir)

    pdata = pData(bg_data)

    mat = log2(as.matrix(gene_fpkms[,-1]) + 1)
    colnames(mat) = gsub('FPKM.','', colnames(mat))

    boxplot_title = 'FPKMs'
    boxplot_y_lab = 'log2(FPKM)'
} else {
    method = 'deseq2'

    results_dir = sprintf('%s/deseq2/02-deseq2_diffex', diffex_dir)
    plots_dir = sprintf('%s/plots', results_dir)

    pdata = data.frame(colData(rld))
    colnames(pdata)[1] = 'sample'

    mat = as.matrix(assay(rld))

    boxplot_title = 'Rlog normalized counts'
    boxplot_y_lab = 'log2(counts)'
}

########################################################
# Plots

########################################################
# PCA and MDS plots for each of the phenotypes

for(phenotype in phenotypes) {

    message(sprintf('Plotting boxplots for %s', phenotype))
    log2_boxplot = plot_boxplot(mat = mat, pdata = pdata, factor_name = phenotype, title = boxplot_title, y_label = boxplot_y_lab, out_name = 'BoxPlot.pdf')

    message(sprintf('Plotting sample heatmap for %s', phenotype))
    log2_heatmap = plot_sample_correlation_heatmap(mat = mat, pdata = pdata, factor_name = phenotype, out_name = 'SampleHeatmap.pdf')

    message(sprintf('Plotting top variably expressed genes heatmap for %s', phenotype))
    plot_top_variably_expressed_heatmap(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 500, out_name = 'Heatmap_TopVar.pdf')

    message(sprintf('Plotting top expressed genes heatmap for %s', phenotype))
    plot_top_expressed_heatmap(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 500, out_name = 'Heatmap_TopExp.pdf')

    # PCA top 500
    message(sprintf('Plotting PCA for %s in dim 1 and 2, top 500', phenotype))
    pca_result_12 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 500, dims = c('PC1','PC2'))
    log2_pca_12 = plot_PCA(compute_PCA_result = pca_result_12, out_name = 'PCAplot_12_top500.pdf')

    message(sprintf('Plotting PCA for %s in dim 2 and 3, top 500', phenotype))
    pca_result_23 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 500, dims = c('PC2','PC3'))
    log2_pca_23 = plot_PCA(compute_PCA_result = pca_result_23, out_name = 'PCAplot_23_top500.pdf')

    message('Plotting scree, top 500')
    scree_plot = plot_scree(compute_PCA_result = pca_result_12, out_name = 'ScreePlot_top500.pdf')

    # PCA top 100
    message(sprintf('Plotting PCA for %s in dim 1 and 2, top 100', phenotype))
    pca_result_12 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 100, dims = c('PC1','PC2'))
    log2_pca_12 = plot_PCA(compute_PCA_result = pca_result_12, out_name = 'PCAplot_12_top100.pdf')

    message(sprintf('Plotting PCA for %s in dim 2 and 3, top 100', phenotype))
    pca_result_23 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 100, dims = c('PC2','PC3'))
    log2_pca_23 = plot_PCA(compute_PCA_result = pca_result_23, out_name = 'PCAplot_23_top100.pdf')

    message('Plotting scree, top 100')
    scree_plot = plot_scree(compute_PCA_result = pca_result_12, out_name = 'ScreePlot_top100.pdf')

    # MDS top 500
    message(sprintf('Plotting MDS for %s in dim 1 and 2, top 500', phenotype))
    log2_mds = plot_MDS(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 500, dims = c(1,2), out_name = 'MDSplot_12_top500.pdf')

    message(sprintf('Plotting MDS for %s in dim 1 and 2, top 500', phenotype))
    log2_mds = plot_MDS(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 500, dims = c(2,3), out_name = 'MDSplot_23_top500.pdf')

    # MDS top 100
    message(sprintf('Plotting MDS for %s in dim 1 and 2, top 500', phenotype))
    log2_mds = plot_MDS(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 100, dims = c(1,2), out_name = 'MDSplot_12_top100.pdf')

    message(sprintf('Plotting MDS for %s in dim 1 and 2, top 500', phenotype))
    log2_mds = plot_MDS(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 100, dims = c(2,3), out_name = 'MDSplot_23_top100.pdf')

}

########################################################
# Volcano plot for each comparison

for(comparison_type in names(de_results)) {

    comparisons = names(de_results[[comparison_type]])

    # Do the gene-level and isoform-level tests for DE
    for(comparison in comparisons) {

        # Extract experiment name and control name for convenience
        out_name = comparison
        exp = unlist(str_split(out_name, '_v_'))[1]
        con = unlist(str_split(out_name, '_v_'))[2]

        message(sprintf('Plotting volcano plot for %s %s', comparison_type, out_name))

        gene_list = data.frame(de_results[[comparison_type]][[comparison]][['gene_list']])

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
