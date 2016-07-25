## abhasi, cgates
## 6/2016
## Watermelon: Recreate Legacy pipeline in snakemake

from collections import defaultdict
import csv
import os

configfile: "config.yaml"

rule all:
    input:
        expand("02-fastqc_reads/{sample}_R1_fastqc.html",
                sample=config["samples"]),
        expand("04-fastqc_align/{sample}_accepted_hits_fastqc.html",
                sample=config["samples"]),
         "05-qc_metrics/Align_summary_all.txt",
#        "06-htseq/HTSeq_counts.txt",
        expand("08-diff_expression/{comparison}/{comparison}_gene.foldchange.{fold_change}.txt",
                comparison=config["comparisons"],
                fold_change=config["fold_change"]),
        expand("08-diff_expression/{comparison}/{comparison}_isoform.foldchange.{fold_change}.txt",
                comparison=config["comparisons"],
                fold_change=config["fold_change"]),
        # expand("09-group_replicates/{comparison}/group_replicates.txt",
        #         comparison=config["comparisons"])
        expand("10-cummerbund/{comparison}/Plots",
                comparison=config["comparisons"]),
        expand("11-deseq/{comparison}/DESeq2_{comparison}_DE.txt",
                comparison=config["comparisons"]),
        expand("11-deseq/{comparison}/DESeq2_{comparison}_DE.txt",
                comparison=config["comparisons"])
        # expand("08-gene_annotation/{comparison}/{comparison}_DEG_ids.txt",
        #          comparison=config["comparisons"]),
        # expand("08-gene_annotation/{comparison}/{comparison}_DEG_names.txt",
        #          comparison=config["comparisons"])
#        expand("10-functional_annotation/{comparison}/{comparison}_DEG_DAVIDchartReport.txt",
#                  comparison=config["comparisons"])

rule concat_reads:
    input:
        "00-multiplexed_reads/{sample}/"
    output:
        "01-raw_reads/{sample}_R1.fastq.gz"
    shell:
        "cat {input}/* > {output}"

rule fastqc_reads:
    input:
        "01-raw_reads/{sample}_R1.fastq.gz"
    output:
        "02-fastqc_reads/{sample}_R1_fastqc.html"
    shell:
        " module load rnaseq && "
        "fastqc {input} -o 02-fastqc_reads"

rule align:
    input:
        "01-raw_reads/{sample}_R1.fastq.gz"
    output:
        "03-tophat/{sample}/{sample}_accepted_hits.bam",
        "03-tophat/{sample}/{sample}_align_summary.txt"
    params:
        gtf_file = config["gtf"],
        transcriptome_index= config["transcriptome_index"],
        bowtie2_index = config["bowtie2_index"],
        sample = lambda wildcards: wildcards.sample
    shell: " module load rnaseq && "
            "tophat -p 8 "
            " --b2-very-sensitive "
            " --no-coverage-search "
            " --library-type fr-unstranded "
            " -I 500000 "
            " -G {params.gtf_file} "
            " --transcriptome-index={params.transcriptome_index} "
            " -T "
            " --no-novel-juncs "
            " -o 03-tophat/{params.sample} "
            " {params.bowtie2_index} "
            " {input} && "
            " mv 03-tophat/{params.sample}/accepted_hits.bam 03-tophat/{params.sample}/{params.sample}_accepted_hits.bam && "
            " mv 03-tophat/{params.sample}/align_summary.txt 03-tophat/{params.sample}/{params.sample}_align_summary.txt "

rule fastqc_align:
    input:
        "03-tophat/{sample}/{sample}_accepted_hits.bam"
    output:
        "04-fastqc_align/{sample}_accepted_hits_fastqc.html"
    shell:
        " module load rnaseq && "
        "fastqc {input} -o 04-fastqc_align"


rule align_qc_metrics:      #Note: this is not the ideal way to do this step; requires "sample_list.txt"; should be simpler
    input:
        expand("03-tophat/{sample}/{sample}_align_summary.txt", sample = config["samples"])
    output:
        "05-qc_metrics/Align_summary_all.txt"
    params:
        align_dir = "03-tophat",
        sample_file = "sample_list.txt",
        output_dir = "05-qc_metrics"
    shell:
        #maybe create sample_list.txt here and then delete it later.
        "perl scripts/getQCMetrics.pl {params.align_dir} {params.sample_file} {params.output_dir}"

rule htseq_per_sample:
    input:
        "03-tophat/{sample}/{sample}_accepted_hits.bam"
    output:
        "06-htseq/{sample}_counts.txt",
    params:
        gtf = config["gtf"],
        input_dir = "06-htseq",
    shell:
        " module load rnaseq &&"
        " python -m HTSeq.scripts.count "
        " -f bam "
        " -s yes "
        " -m "
        " intersection-nonempty "
        " -q {input} "
        " {params.gtf} "
        " > {output} && "
        "perl scripts/mergeHTSeqCountFiles.pl {params.input_dir}"

# unable to do this in one step; gives error: 'Not all output files of rule htseq_per_sample contain the same wildcards.'
rule htseq_merge:
    input:
        expand("06-htseq/{sample}_counts.txt", sample=config["samples"])
    output:
        "06-htseq/HTSeq_counts.txt"
    params:
        output_dir = "06-htseq",
        input_dir = "06-htseq"
    shell:
       " perl scripts/mergeHTSeqCountFiles.pl {params.input_dir} "

def cuffdiff_labels(underbar_separated_comparisons):
    return underbar_separated_comparisons.replace("_", ",")

def cuffdiff_samples(underbar_separated_comparisons,
                     sample_name_group,
                     sample_file_format):
    group_sample_names = defaultdict(list)
    for actual_sample_name, group in sample_name_group.items():
        group_sample_names[group].append(sample_file_format.format(sample_placeholder=actual_sample_name))

    params = []
    for group in underbar_separated_comparisons.split('_'):
        params.append(','.join(sorted(group_sample_names[group])))

    return ' '.join(params)

rule cuffdiff:
    input:
        expand("03-tophat/{sample}/{sample}_accepted_hits.bam", sample=config["samples"])
    output:
        "07-cuffdiff/{comparison}/gene_exp.diff",
        "07-cuffdiff/{comparison}/isoform_exp.diff",
        "07-cuffdiff/{comparison}/read_groups.info"
    params:
        output_dir = "07-cuffdiff/{comparison}",
        labels = lambda wildcards : cuffdiff_labels(wildcards.comparison),
        samples = lambda wildcards : cuffdiff_samples(wildcards.comparison,
                                                      config["samples"],
                                                      "03-tophat/{sample_placeholder}/{sample_placeholder}_accepted_hits.bam"),
        fasta = config["fasta"],
        gtf_file = config["gtf"],
    shell:
        " module load rnaseq && "
        " cuffdiff -q "
        " -L {params.labels} "  # Wt, Mut1, Mut2
        " --max-bundle-frags 999999999 "
        " -o {params.output_dir} "
        " -b {params.fasta} "
        " -u -N --compatible-hits-norm "
        " {params.gtf_file} "
        " {params.samples}" #Sample_50153_accepted_hits.bam Sample_50154_accepted_hits.bam Sample_50155_accepted_hits.bam .....

rule diff_exp:
    input:
        cuffdiff_gene_exp="07-cuffdiff/{comparison}/gene_exp.diff",
        cuffdiff_isoform_exp="07-cuffdiff/{comparison}/isoform_exp.diff"
    output:
        "08-diff_expression/{comparison}/{comparison}_gene.foldchange.{fold_change}.txt",
        "08-diff_expression/{comparison}/{comparison}_isoform.foldchange.{fold_change}.txt"
    params:
        output_dir = "08-diff_expression/{comparison}",
        fold_change = config["fold_change"],
        comparison = lambda wildcards: wildcards.comparison
    shell:
        " module load rnaseq && "
        " mkdir -p {params.output_dir} &&"
        " perl scripts/Cuffdiff_out_format_v6.pl "
        " {input.cuffdiff_gene_exp} "
        " {params.comparison}_gene "
        " {params.fold_change} "
        " {params.output_dir} && "
        " perl scripts/Cuffdiff_out_format_v6.pl "
        " {input.cuffdiff_isoform_exp} "
        " {params.comparison}_isoform "
        " {params.fold_change} "
        " {params.output_dir} "

rule deseq_build_group_replicates:
    input:
        "07-cuffdiff/{comparison}/read_groups.info"
    output:
        "09-group_replicates/{comparison}/group_replicates.txt"
    params:
        comparison = lambda wildcards: wildcards.comparison
    run:
        input_file_name = input[0]
        output_file_name = output[0]
        replicate_dict = defaultdict(list)

        with open(input_file_name, "r") as input_file, \
            open(output_file_name, "w") as output_file:
            reader = csv.DictReader(input_file, delimiter='\t')
            for row in reader:
                sample = os.path.basename(row["file"]).split("_accepted_hits.bam")[0]
                replicate_dict[row["condition"]].append("replicate_{}: {}".format(row["replicate_num"], sample))

            for group, samples in replicate_dict.items():
                all_samples = ", ".join(samples)
                output_file.write("{}\t{}\n".format(group, all_samples))

rule cummerbund:
    input:
        "09-group_replicates/{comparison}/group_replicates.txt"
    output:
        "10-cummerbund/{comparison}/Plots"
    params:
        cuff_diff_dir = "07-cuffdiff/{comparison}",
        output_dir = "10-cummerbund/{comparison}/",
        gtf_file = config["gtf"],
        genome = config["genome"],
#        logfile = "10-cummerbund/{comparison}/cummerbund.log"
    log:
         "10-cummerbund/{comparison}/cummerbund.log"
    shell:
        " module load rnaseq && "
        " mkdir -p {output} && "
        " Rscript scripts/Run_cummeRbund_v1.R "
        " baseDir={params.output_dir} "
        " cuffDiffDir=../../{params.cuff_diff_dir} "
        " grpRepFile={input} "
        " gtfFile={params.gtf_file} "
        " genome={params.genome} --log {log}"

rule deseq2:
    input:
        counts_file = "06-htseq/HTSeq_counts.txt",
        group_replicates = "09-group_replicates/{comparison}/group_replicates.txt"
    output:
        "11-deseq/{comparison}/DESeq2_{comparison}_DE.txt",
        "11-deseq/{comparison}/DESeq2_{comparison}_DESig.txt"
    params:
        output_dir = "11-deseq",
        comparison = lambda wildcards: wildcards.comparison
    shell:
        " module load rnaseq && "
        " Rscript scripts/Run_DESeq.R "
        " baseDir={params.output_dir}/{params.comparison} "
        " grpRepFile={input.group_replicates} "
        " countsFile={input.counts_file} "

# rule gene_annotation:
#     input:
#         "07-cuffdiff/{comparison}/isoform_exp.diff"
#     output:
#         ids = "08-gene_annotation/{comparison}_DEG_ids.txt",
#         names = "08-gene_annotation/{comparison}_DEG_names.txt"
#     params:
#         output_dir = "08-gene_annotation",
#         input_dir = "07-cuffdiff",
#         genome = "hg19",
#         gtf_file = config["gtf"],
#         info_file = "info.txt",
#         comparison = lambda wildcards: wildcards.comparison
#     shell:
#         "touch {params.output_dir}/{params.info_file} && "
#         "perl scripts/Jess_scripts/get_NCBI_gene_annotation_JESS.pl {params.input_dir}/{params.comparison} {params.output_dir} {params.genome} {params.gtf_file} {params.output_dir} {params.output_dir}/{params.info_file} && "
#         "cp {params.input_dir}/{params.comparison}/{params.comparison}_DEG* {params.output_dir} && "
#         "touch {output.ids} {output.names}" #have to touch output because snakemake gets confused about the timestamps when an archive is extracted
#
# rule functional_annotation:
#     input:
#         "09-gene_annotation/{comparison}/{comparison}_DEG_ids.txt"
#     output:
#         report = "10-functional_annotation/{comparison}/{comparison}_DEG_DAVIDchartReport.txt"
#     params:
#         input_dir = "09-gene_annotation/{comparison}",
#         output_dir = "10-functional_annotation",
#         comparison = lambda wildcards: wildcards.comparison
#     shell:
#         " module load rnaseq && "
#         "perl scripts/get_DAVIDchartReport.pl -dir {params.input_dir} -type ENTREZ_GENE_ID -thd 0.05 && "
#         "cp {params.input_dir}/{params.comparison}_DEG_DAVIDchartReport.txt {params.output_dir} && "
#         "touch {output} "
