library(optparse)

option_list = list(
    make_option('--stringtie_dir', type='character', help='[Required] Path to working directory'),
    make_option('--config_file', type='character', help='[Required] Path to configfile')
)
opt = parse_args(OptionParser(option_list=option_list))

stringtie_dir = opt$stringtie_dir
config_file = opt$config_file

#######################################

library(yaml)
library(stringr)
library(ballgown)

#######################################

carat_cleanup = function(str) {
    tmp = unlist(strsplit(str, '^', fixed = TRUE))
    tmp = str_trim(tmp)
    return(tmp)
}

versus_split = function(str) {
    tmp = unlist(strsplit(str, '_v_', fixed = TRUE))
    tmp = str_trim(tmp)
    return(tmp)
}

handle_phenotype_data = function(yaml, sample_paths) {
    yaml_phenotypes = carat_cleanup(yaml$phenotypes)

    yaml_samples = names(yaml$samples)

    yaml_phenotype_data = sapply(yaml_samples,
        function(sample){
            carat_cleanup(yaml$samples[[sample]])
        }, USE.NAMES = FALSE)

    # If there is only one phenotype attribute per sample
    if(class(yaml_phenotype_data) == 'character') {
        df = data.frame(yaml_phenotype_data)
        colnames(df) = yaml_phenotypes
    } else {
        df = data.frame(t(yaml_phenotype_data))
        colnames(df) = yaml_phenotypes
    }

    # Construct the pdata object by adding the column of sample names and match
    # them in the order of the sample_paths
    pdata = data.frame(
        sample = yaml_samples,
        df
    )
    pdata = pdata[match(basename(sample_paths), pdata$sample), ]

    return(pdata)
}

handle_comparisons = function(yaml) {
    comps = lapply(yaml$comparisons,
        function(comparison){
            tmp = data.frame(t(sapply(comparison, versus_split, USE.NAMES = FALSE)))
            colnames(tmp) = c('exp', 'con')
            return(tmp)
        }
    )
    return(comps)
}

#######################################
# yaml parsing
yaml = read_yaml(config_file)

diffex_dir = yaml$diffex_output_dir
results_dir = sprintf('%s/ballgown/01-ballgown_diffex', diffex_dir)
ballgown_inputs_dir = sprintf('%s/ballgown', stringtie_dir)

sample_paths = list.files(ballgown_inputs_dir, full.names = TRUE)

if(length(sample_paths) == 0) {
    stop(sprintf('No directories listed in %s', ballgown_inputs_dir))
}

pdata = handle_phenotype_data(yaml, sample_paths)

comparisons = handle_comparisons(yaml)

comparison_types = names(comparisons)

fdr_cutoff = yaml$deseq2_adjustedPValue
fc_cutoff = yaml$fold_change

#######################################
message('Reading ballgown input')
bg_data = ballgown(samples = sample_paths, meas = 'all', pData = pdata)

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

#######################################
# Loop through the comparison types, and the comparisons
for(type in comparison_types) {
    message(sprintf('On %s comparison type', type))
    # What are all the comparisons to be done based on the type?
    type_comparisons = comparisons[[type]]

    # Setup the results directories
    de_dir = sprintf('%s/%s', results_dir, type)

    # Do the gene-level and isoform-level tests for DE
    for(i in 1:nrow(type_comparisons)) {
        # Extract experiment name and control name for convenience
        exp = as.character(type_comparisons[i, 'exp'])
        con = as.character(type_comparisons[i, 'con'])
        out_name = paste(exp, con, sep='_v_')

        message(sprintf('On %s vs %s comparison', exp, con))
        # Wow this is BAD. Why not just extend SummarizedExperiment?
        sub_bg_data = subset(bg_data, sprintf('%s == "%s" | %s == "%s"', type, exp, type, con), genomesubset = FALSE)
        # Relevel the factor so the reference is correct
        pData(sub_bg_data)[, type] = relevel(pData(sub_bg_data)[, type], ref = con)

        gene_results = stattest(
            gown = sub_bg_data,
            feature = 'gene',
            meas='FPKM',
            covariate = type,
            getFC = TRUE)
        # Order by qval, and remove the feature column
        gene_results = gene_results[order(gene_results$qval, gene_results$fc), c('id', 'fc', 'pval', 'qval')]
        gene_results$Condition = exp
        gene_results$Control = con
        gene_results$diff_exp = ifelse(gene_results$fc > fc_cutoff & gene_results$qval < fdr_cutoff, 'Yes', 'No')

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
        iso_results$Condition = exp
        iso_results$Control = con
        iso_results$diff_exp = ifelse(iso_results$fc > fc_cutoff & iso_results$qval < fdr_cutoff, 'Yes', 'No')

        write.table(gene_results, file = sprintf('%s/%s_gene.txt', de_dir, out_name), sep='\t', col.names = T, row.names = F, quote = F)
        write.table(iso_results, file = sprintf('%s/%s_isoform.txt', de_dir, out_name), sep='\t', col.names = T, row.names = F, quote = F)
    }
}

# message('Saving ballgown RData')
# save(bg_data, file = sprintf('%s/ballgown_data.rda', results_dir), compress = 'xz')
