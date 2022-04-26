##########
# Set up logging and save snakemake S4 object (for debugging or running manually)
log = file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')
save(snakemake, file = snakemake@params[['snakemake_rdata']])

#setwd("/nfs/med-bfx-activeprojects/trsaari/sandbox/20191023_testing_simulated_data/")
#load("analysis_20191024/diffex_results/plots_labeled_by_pheno/.deseq2_plots_by_phenotype_snakemake.rda")

#Isolate conda environment: https://github.com/conda-forge/r-base-feedstock/issues/37
#If we move away from conda in the future, we may want to remove this
.libPaths(R.home("library"))

########################################################
# Load libraries
lib.vector = c("DESeq2", "dplyr", "GGally", "ggfortify", "ggplot2",
    "ggrepel", "pheatmap", "RColorBrewer", "stringr", "tidyr", "yaml")

foo = suppressMessages(lapply(lib.vector, library, character.only=T, warn.conflicts=F, quietly=T))

########################################################
# Define plotting functions

make_heatmap_annots = function(pdata, factors) {
    if(!is.null(factors)) {
        if(!all(factors %in% colnames(pdata))) {
            stop('Not all factor_covs are in colnames(pdata).')
        }

        annot_df = data.frame(pdata[, factors])
        # If there is only one factor, we have to coerce to data.frame
        # but there is no harm when more than one factor. Moreover,
        # we have to supply the correct column and row names in this case
        if(ncol(annot_df) == 1) {
            rownames(annot_df) = rownames(pdata)
            colnames(annot_df) = factors
        }

        annot_colors = lapply(factors, function(factor){
            factor_levels = sort(unique(pdata[, factor]))
            if(length(factor_levels) <= 9){
                color_vector = RColorBrewer::brewer.pal(length(factor_levels), 'Set1')
            } else if(length(factor_levels) > 9 && length(factor_levels) <= 13){
                color_vector = RColorBrewer::brewer.pal(length(factor_levels), 'Set3')
            }

            if(length(factor_levels) > length(color_vector)){
                msg = paste0('Too many factor levels in column ', factor, ', cannot assign color palette')
                stop(msg)
            } else if(length(factor_levels) < length(color_vector)) {
                # If there are more colors than needed (minimum palette size is 3)
                # then shrink color vector to appropiate size
                color_vector = color_vector[1:length(factor_levels)]
            }
            names(color_vector) = factor_levels

            return(color_vector)
        })
        names(annot_colors) = factors
    } else {
        annot_df = NA
        annot_colors = NA
    }

    return(list(annot_df = annot_df, annot_colors = annot_colors))
}


plot_sample_correlation_heatmap = function(mat, pdata, factors, out_basename = 'SampleHeatmap') {

    annot_list = make_heatmap_annots(pdata, factors)

    # Calculate the distance matrix and clustering object
    dist_obj = dist(t(mat), method = 'euclidean', diag = T, upper = T)
    hclust_obj = hclust(dist_obj, method = 'complete')
    dist_mat = as.matrix(dist_obj)

    pheatmap( # pdf
        mat = dist_mat,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        angle_col = 45,
        display_numbers = TRUE,
        cellwidth = 20,
        cellheight = 20,
        fontsize = 8,
        fontsize_row = 6,
        fontsize_col = 6,
        fontsize_number = 6,
        number_color = 'black',
        annotation_col = annot_list[['annot_df']],
        annotation_colors = annot_list[['annot_colors']],
        color = colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))(20),
        filename = file.path(plots_dir, paste0(out_basename, '.pdf'))
    )
    pheatmap( # png
        mat = dist_mat,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        angle_col = 45,
        display_numbers = TRUE,
        cellwidth = 20,
        cellheight = 20,
        fontsize = 8,
        fontsize_row = 6,
        fontsize_col = 6,
        fontsize_number = 6,
        number_color = 'black',
        annotation_col = annot_list[['annot_df']],
        annotation_colors = annot_list[['annot_colors']],
        color = colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))(20),
        filename = file.path(plots_dir, paste0(out_basename, '.png'))
    )

    return(list(dist_obj = dist_obj, hclust_obj = hclust_obj))
}


plot_top_variably_expressed_heatmap = function(mat, pdata, factors, top_n = 1000, out_basename = 'Heatmap_TopVar') {

    annot_list = make_heatmap_annots(pdata, factors)

    #Add guard for when mat has less rows than specified top_n
    top_n = min(nrow(mat), top_n)

    # Calculate the top_n variable genes and put them in decreasing order
    top_var_mat = mat[order(matrixStats::rowVars(mat), decreasing = T), ][1:top_n, ]

    pdf(file = file.path(plots_dir, paste0(out_basename, '.pdf')), height = 20, width = 10)
        pheatmap(
            mat = top_var_mat,
            scale= 'row',
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = FALSE,
            annotation_col = annot_list[['annot_df']],
            annotation_colors = annot_list[['annot_colors']],
            fontsize = 7,
            fontsize_row = 7,
            las = 2,
            main = sprintf('Top %s variably expressed genes', top_n),
            color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Blues')))(255)
        )
    dev.off()
    png(file = file.path(plots_dir, paste0(out_basename, '.png')), height = 2400, width = 1200)
        pheatmap(
            mat = top_var_mat,
            scale= 'row',
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = FALSE,
            annotation_col = annot_list[['annot_df']],
            annotation_colors = annot_list[['annot_colors']],
            fontsize = 7,
            fontsize_row = 7,
            las = 2,
            main = sprintf('Top %s variably expressed genes', top_n),
            color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Blues')))(255)
        )
    dev.off()
}


plot_top_expressed_heatmap = function(mat, pdata, factors, top_n = 1000, out_basename = 'Heatmap_TopExp') {

    annot_list = make_heatmap_annots(pdata, factors)

    #Add guard for when mat has less rows than specified top_n
    top_n = min(nrow(mat), top_n)

    # Calculate the top_n expressed genes and put them in decreasing order
    top_exp_mat = mat[order(rowMeans(mat), decreasing = T), ][1:top_n, ]

    pdf(file = file.path(plots_dir, paste0(out_basename, '.pdf')), height = 20, width = 10)
        pheatmap(
            mat = top_exp_mat,
            scale= 'row',
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = FALSE,
            annotation_col = annot_list[['annot_df']],
            annotation_colors = annot_list[['annot_colors']],
            fontsize = 7,
            fontsize_row = 7,
            las = 2,
            main = sprintf('Top %s expressed genes', top_n),
            color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Blues')))(255)
        )
    dev.off()
    png(file = file.path(plots_dir, paste0(out_basename, '.png')), height = 2400, width = 1200)
        pheatmap(
            mat = top_exp_mat,
            scale= 'row',
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = FALSE,
            annotation_col = annot_list[['annot_df']],
            annotation_colors = annot_list[['annot_colors']],
            fontsize = 7,
            fontsize_row = 7,
            las = 2,
            main = sprintf('Top %s expressed genes', top_n),
            color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Blues')))(255)
        )
    dev.off()
}

plot_boxplot = function(mat, pdata, factor_name, title, y_label, out_basename = 'BoxPlot') {

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
    ggsave(filename = file.path(plots_dir, factor_name, paste0(out_basename, '.pdf')), plot = box_plot, height = 8, width = 8, dpi = 300)
    ggsave(filename = file.path(plots_dir, factor_name, paste0(out_basename, '.png')), plot = box_plot, height = 8, width = 8, dpi = 300)

    return(box_plot)
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

    #Add guard for when mat has less rows than specified top_n
    top_n = min(nrow(mat), top_n)

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


plot_scree = function(compute_PCA_result, out_basename = 'ScreePlot') {
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
    ggsave(filename = file.path(plots_dir, factor_name, paste0(out_basename, '.pdf')), plot = scree_plot, height = 6, width = 6, dpi = 300)
    ggsave(filename = file.path(plots_dir, factor_name, paste0(out_basename, '.png')), plot = scree_plot, height = 6, width = 6, dpi = 300)

    return(scree_plot)
}


plot_PCA = function(compute_PCA_result, out_basename = 'PCAplot') {

    pca_df = compute_PCA_result[['pca_df']]
    var_explained = compute_PCA_result[['var_explained']]
    top_n = compute_PCA_result[['top_n']]
    factor_name = compute_PCA_result[['factor_name']]
    groups = as.factor(pca_df[,factor_name])
    replicates = as.factor(pca_df[,'replicate'])
    dims = compute_PCA_result[['dims']]

    # Plot PCA
    # Use colors for groups or reps, whichever has more numerous levels
    if(nlevels(groups) > 6 && nlevels(replicates) > 6) {
      message(paste0("Warning - not enough symbols to represent all groups and replicates. Assigning colors to groups and ignoring replicates"))
      pca_plot = ggplot(pca_df, aes(x = x, y = y, color = groups)) +
        geom_point(alpha=0.6) +
        labs(
          title = sprintf('%s PCA plot', factor_name),
          subtitle = sprintf('Using top %s variable genes', top_n),
          x = sprintf('%s: %s%% variance', dims[1], var_explained[1]),
          y = sprintf('%s: %s%% variance', dims[2], var_explained[2])) +
        theme_bw()
    } else if(nlevels(groups) > nlevels(replicates)){
      pca_plot = ggplot(pca_df, aes(x = x, y = y, color = groups)) +
        geom_point(aes_string(shape = 'replicate')) +
        labs(
          title = sprintf('%s PCA plot', factor_name),
          subtitle = sprintf('Using top %s variable genes', top_n),
          x = sprintf('%s: %s%% variance', dims[1], var_explained[1]),
          y = sprintf('%s: %s%% variance', dims[2], var_explained[2])) +
        theme_bw()
    } else {
      pca_plot = ggplot(pca_df, aes(x = x, y = y, color = replicates)) +
        geom_point(aes_string(shape = factor_name)) +
        labs(
          title = sprintf('%s PCA plot', factor_name),
          subtitle = sprintf('Using top %s variable genes', top_n),
          x = sprintf('%s: %s%% variance', dims[1], var_explained[1]),
          y = sprintf('%s: %s%% variance', dims[2], var_explained[2])) +
        theme_bw()
    }

    ggsave(filename = file.path(plots_dir, factor_name, paste0(out_basename, '.pdf')), plot = pca_plot, height = 6, width = 6, dpi = 300)
    ggsave(filename = file.path(plots_dir, factor_name, paste0(out_basename, '.png')), plot = pca_plot, height = 6, width = 6, dpi = 300)

    return(pca_plot)
}


########################################################
########################################################
########################################################

phenotypes = snakemake@params[['phenotypes']]

diffex_dir = snakemake@params[['diffex_dir']]
#Remove trailing / if there is one
diffex_dir = sub("/$", "", diffex_dir)

########################################################
# Load data

load(snakemake@input[['rda']])

if('bg_data' %in% ls()) {
    method = 'ballgown'

    plots_dir = sprintf('%s/ballgown/plots', diffex_dir)

    pdata = pData(bg_data)

    mat = log2(as.matrix(gene_fpkms[,-1]) + 1)
    colnames(mat) = gsub('FPKM.','', colnames(mat))

    boxplot_title = 'FPKMs'
    boxplot_y_lab = 'log2(FPKM)'
} else {
    method = 'deseq2'

    plots_dir = sprintf('%s/plots_labeled_by_pheno', diffex_dir)

    pdata = data.frame(colData(rld))
    colnames(pdata)[1] = 'sample'
    #These are the strings we want. If sample names are 001, 002, etc, they are carried forward this way
    pdata$sample = rownames(pdata)

    mat = as.matrix(assay(rld))

    boxplot_title = 'Rlog normalized counts'
    boxplot_y_lab = 'log2(counts)'
}

########################################################
# Plots

########################################################
# Plots with all phenotypes included

message(sprintf('Plotting sample heatmap for %s', paste(phenotypes, collapse = ', ')))
log2_heatmap = plot_sample_correlation_heatmap(mat = mat, pdata = pdata, factors = phenotypes, out_basename = 'SampleHeatmap')

message(sprintf('Plotting top variably expressed genes heatmap for %s', paste(phenotypes, collapse = ', ')))
plot_top_variably_expressed_heatmap(mat = mat, pdata = pdata, factors = phenotypes, top_n = 500, out_basename = 'Heatmap_TopVar')

message(sprintf('Plotting top expressed genes heatmap for %s', paste(phenotypes, collapse = ', ')))
plot_top_expressed_heatmap(mat = mat, pdata = pdata, factors = phenotypes, top_n = 500, out_basename = 'Heatmap_TopExp')

# Plots for each of the phenotypes
for(phenotype in phenotypes) {

    message(sprintf('Plotting boxplots for %s', phenotype))
    log2_boxplot = plot_boxplot(mat = mat, pdata = pdata, factor_name = phenotype, title = boxplot_title, y_label = boxplot_y_lab, out_basename = 'BoxPlot_rlog')
    raw_boxplot = plot_boxplot(mat = log2(raw_counts), pdata = pdata, factor_name = phenotype, title = 'Non-normalized counts', y_label = boxplot_y_lab, out_basename = 'BoxPlot_raw')

    # PCA top 500
    message(sprintf('Plotting PCA for %s in dim 1 and 2, top 500', phenotype))
    pca_result_12 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 500, dims = c('PC1','PC2'))
    log2_pca_12 = plot_PCA(compute_PCA_result = pca_result_12, out_basename = 'PCAplot_12_top500')

    message(sprintf('Plotting PCA for %s in dim 2 and 3, top 500', phenotype))
    pca_result_23 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 500, dims = c('PC2','PC3'))
    log2_pca_23 = plot_PCA(compute_PCA_result = pca_result_23, out_basename = 'PCAplot_23_top500')

    message('Plotting scree, top 500')
    scree_plot = plot_scree(compute_PCA_result = pca_result_12, out_basename = 'ScreePlot_top500')

    # PCA top 100
    message(sprintf('Plotting PCA for %s in dim 1 and 2, top 100', phenotype))
    pca_result_12 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 100, dims = c('PC1','PC2'))
    log2_pca_12 = plot_PCA(compute_PCA_result = pca_result_12, out_basename = 'PCAplot_12_top100')

    message(sprintf('Plotting PCA for %s in dim 2 and 3, top 100', phenotype))
    pca_result_23 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 100, dims = c('PC2','PC3'))
    log2_pca_23 = plot_PCA(compute_PCA_result = pca_result_23, out_basename = 'PCAplot_23_top100')

    message('Plotting scree, top 100')
    scree_plot = plot_scree(compute_PCA_result = pca_result_12, out_basename = 'ScreePlot_top100')

}
