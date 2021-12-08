##########
# Set up logging and save snakemake S4 object (for debugging or running manually)
log = file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')
save(snakemake, file = snakemake@params[['snakemake_rdata']])

#setwd("/nfs/med-bfx-activeprojects/trsaari/sandbox/20191023_testing_simulated_data/")
#load("analysis_20191024/diffex_results/plots_labeled_by_pheno/.deseq2_heatmaps_snakemake.rda")

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
            color_vector = RColorBrewer::brewer.pal(length(factor_levels), 'Set1')
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


########################################################
########################################################
########################################################

phenotypes = names(snakemake@params[['sample_phenotypes']])

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
} else {
    method = 'deseq2'

    plots_dir = sprintf('%s/plots_labeled_by_pheno', diffex_dir)

    pdata = data.frame(colData(rld))
    colnames(pdata)[1] = 'sample'
    #These are the strings we want. If sample names are 001, 002, etc, they are carried forward this way
    pdata$sample = rownames(pdata)

    mat = as.matrix(assay(rld))
}

########################################################
# Plots

message(sprintf('Plotting sample heatmap for %s', paste(phenotypes, collapse = ', ')))
log2_heatmap = plot_sample_correlation_heatmap(mat = mat, pdata = pdata, factors = phenotypes, out_basename = 'SampleHeatmap')

message(sprintf('Plotting top variably expressed genes heatmap for %s', paste(phenotypes, collapse = ', ')))
plot_top_variably_expressed_heatmap(mat = mat, pdata = pdata, factors = phenotypes, top_n = 500, out_basename = 'Heatmap_TopVar')

message(sprintf('Plotting top expressed genes heatmap for %s', paste(phenotypes, collapse = ', ')))
plot_top_expressed_heatmap(mat = mat, pdata = pdata, factors = phenotypes, top_n = 500, out_basename = 'Heatmap_TopExp')
