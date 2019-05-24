library(optparse)

option_list = list(
    make_option('--rsem_dir', type='character', help='[Required] Path to working directory'),
    make_option('--config_file', type='character', help='[Required] Path to configfile')
)
opt = parse_args(OptionParser(option_list=option_list))

rsem_dir = opt$rsem_dir
#config_file = opt$config_file

rsem_dir = "/nfs/med-bfx-activeprojects/trsaari/example_output_Watermelon/alignment_results/04-rsem_star_align" #TWS DEBUG
config_file = "/nfs/med-bfx-activeprojects/trsaari/Watermelon/config/config_190521_test.yaml" #TWS DEBUG

#######################################

library(ballgown)
library(readr)
library(stringr)
library(yaml)

#######################################

versus_split = function(str) {
    tmp = unlist(strsplit(str, '_v_', fixed = TRUE))
    tmp = str_trim(tmp)
    return(tmp)
}

handle_comparisons = function(yaml) {
    comps = lapply(yaml$comparisons,
        function(comparison){
            tmp = data.frame(t(sapply(comparison, versus_split, USE.NAMES = FALSE)))
            colnames(tmp) = c('exp', 'con')
            tmp$out_name = apply(tmp, 1, paste, collapse = '_v_')
            return(tmp)
        }
    )
    return(comps)
}

get.samplenames.from.dir <- function(results.dir) {
  fnames <- list.files(results.dir)
  iso.fnames <- fnames[grepl(".isoforms.results", fnames)]
  sample.names <- sub(".isoforms.results", "", iso.fnames, fixed=T)
  return(sample.names)
}

#######################################
# yaml parsing

yaml = read_yaml(config_file)
#diffex_yaml = read_yaml(yaml$diffex$comparison_file)
diffex_yaml = read_yaml("/nfs/med-bfx-activeprojects/trsaari/Watermelon/config/example_comparisons.yaml") #TWS DEBUG

# Establish directories
diffex_dir = yaml$dirs$diffex_output
results_dir = sprintf('%s/ballgown/01-ballgown_diffex', diffex_dir)
counts_dir = sprintf('%s/counts', results_dir)
gene_lists_dir = sprintf('%s/gene_lists', results_dir)
#ballgown_inputs_dir = sprintf('%s/ballgown', stringtie_dir)
ballgown_inputs_dir = sprintf('%s/ballgown', rsem_dir)

# Get phenotyp matrix
#pdata = read_csv(yaml$diffex$sample_description_file)
pdata = read_csv("/nfs/med-bfx-activeprojects/trsaari/Watermelon/config/example_sample_description_subset.csv") #TWS DEBUG

# Establish comparisons
comparisons = handle_comparisons(diffex_yaml)
comparison_types = names(comparisons)

# Establish cutoffs
fdr_cutoff = diffex_yaml$deseq2_adjustedPValue
fc_cutoff = log2(diffex_yaml$fold_change)

#######################################
#message('Reading ballgown input')
# sample_paths = list.files(ballgown_inputs_dir, full.names = TRUE) #TWS - This isn't needed for rsem input

# if(length(sample_paths) == 0) {
#     stop(sprintf('No directories listed in %s', ballgown_inputs_dir))
# }

yaml$references$gtf <- "/nfs/med-bfx-activeprojects/trsaari/sandbox/hg19_noRibo.gtf" #TWS DEBUG
pdata = as.data.frame(pdata) #TWS DEBUG

rsem_samplenames <- get.samplenames.from.dir(rsem_dir)
bg_data = ballgownrsem(dir=rsem_dir, samples=rsem_samplenames, gtf=yaml$references$gtf, meas = 'all', pData = pdata)

###################
# Pull out gene-level FPKMs
gene_fpkms = gexpr(bg_data)
# Retain the id column (gene symbols)
gene_fpkms = cbind(id = rownames(gene_fpkms), data.frame(gene_fpkms))

###################
# Pull out transcript-level FPKMs
iso_fpkms = texpr(bg_data, meas = 'all')

# Create a dictionary between t_id (internal to ballgown), gene_id (symbol),
# the t_name (transcript RefSeq name), and transcript location / strand
# NOTE: This information is determined from the gtf files input to ballgown()
# which are in turn determined based on the gtf given in the config file.
gtf_dict = iso_fpkms[, c('t_id', 'gene_id', 't_name', 'chr', 'start', 'end', 'strand')]

# Retain identifying information for the transcripts
iso_col_keep = c('t_id', 'gene_id', 't_name', 'chr', 'start', 'end', 'strand', grep('FPKM', colnames(iso_fpkms), value = T))
iso_fpkms = iso_fpkms[, iso_col_keep]

###################
# Write both to a tab-delimited file
gene_fpkm_file = sprintf('%s/gene_fpkms.txt', counts_dir)
iso_fpkm_file = sprintf('%s/iso_fpkms.txt', counts_dir)

message(sprintf('Writing gene FPKMs to %s', gene_fpkm_file))
write.table(gene_fpkms, file = gene_fpkm_file, sep = '\t', row.names = F, col.names = T, quote = F)

message(sprintf('Writing isoform FPKMs to %s', iso_fpkm_file))
write.table(iso_fpkms, file = iso_fpkm_file, sep = '\t', row.names = F, col.names = T, quote = F)

#######################################
# Loop through the comparison types, and the comparisons

# List container for all the results
de_results = list()

for(type in comparison_types) {
    message(sprintf('On %s comparison type', type))
    # What are all the comparisons to be done based on the type?
    type_comparisons = comparisons[[type]]

    # Setup the results directories
    de_dir = sprintf('%s/%s', gene_lists_dir, type)

    # List container named type_comparisons with gene and isoform sublists
    type_de_results = list()

    # Do the gene-level and isoform-level tests for DE
    for(i in 1:nrow(type_comparisons)) {
        # Extract experiment name and control name for convenience
        exp = as.character(type_comparisons[i, 'exp'])
        con = as.character(type_comparisons[i, 'con'])
        out_name = as.character(type_comparisons[i, 'out_name'])

        type_de_results[[out_name]] = list()

        message(sprintf('On %s vs %s comparison', exp, con))
        # Wow this is BAD. Why not just extend SummarizedExperiment?
        sub_bg_data = subset(bg_data, sprintf('%s == "%s" | %s == "%s"', type, exp, type, con), genomesubset = FALSE)
        # Relevel the factor so the reference is correct
        pData(sub_bg_data)[, type] = relevel(pData(sub_bg_data)[, type], ref = con)
        # NOTE: We have subsetted the bg_data, so the library normalization, etc. is on the subset
        # Because there is no notion of contrast in ballgown, there is no way to

        gene_results = stattest(
            gown = sub_bg_data,
            feature = 'gene',
            meas='FPKM',
            covariate = type,
            getFC = TRUE)
        # Order by qval, and remove the feature column
        gene_results = gene_results[order(gene_results$qval, gene_results$fc), c('id', 'fc', 'pval', 'qval')]
        gene_results$log2fc = log2(gene_results$fc)
        gene_results$Condition = exp
        gene_results$Control = con
        gene_results$diff_exp = ifelse(abs(gene_results$log2fc) >= fc_cutoff & gene_results$qval < fdr_cutoff, 'YES', 'NO')

        iso_results = stattest(
            gown = sub_bg_data,
            feature = 'transcript',
            meas='FPKM',
            covariate = type,
            getFC = TRUE)
        # Annotate the isoform results with the gene_id and t_name from iso_fpkms
        iso_results = merge(
            x = gtf_dict,
            y = iso_results,
            by.x = 't_id',
            by.y = 'id',
            all.x = FALSE,
            all.y = TRUE)
        # Order by qval, and remove the feature and t_id columns
        iso_results = iso_results[order(iso_results$qval, iso_results$fc), setdiff(colnames(iso_results), c('t_id','feature'))]
        iso_results$log2fc = log2(iso_results$fc)
        iso_results$Condition = exp
        iso_results$Control = con

        iso_results$diff_exp = ifelse(abs(iso_results$log2fc) >= fc_cutoff & iso_results$qval < fdr_cutoff, 'YES', 'NO')

        # Sort by diff_exp, then qval
        gene_results = gene_results[order(-rank(gene_results$diff_exp), gene_results$qval), ]
        iso_results = iso_results[order(-rank(iso_results$diff_exp), iso_results$qval), ]

        # Append to the type_de_results
        type_de_results[[out_name]] = list(gene_list = gene_results, isoform_list = iso_results)

        gene_list_file = sprintf('%s/%s_gene.txt', de_dir, out_name)
        message(sprintf('Writing gene list to %s', gene_list_file))
        write.table(gene_results, file = gene_list_file, sep='\t', col.names = T, row.names = F, quote = F)

        iso_list_file = sprintf('%s/%s_isoform.txt', de_dir, out_name)
        message(sprintf('Writing isoform list to %s', iso_list_file))
        write.table(iso_results, file = iso_list_file, sep='\t', col.names = T, row.names = F, quote = F)
    }

    de_results[[type]] = type_de_results
}

# Save the rda file that gets passed on to plotting rule
ballgown_rda_file = sprintf('%s/ballgown_data.rda', results_dir)
save_list = c('yaml', 'comparisons', 'comparison_types', 'bg_data', 'gene_fpkms', 'iso_fpkms', 'gtf_dict', 'de_results')

message(sprintf('Saving ballgown RData to %s', ballgown_rda_file))
save(list = save_list, file = ballgown_rda_file)
