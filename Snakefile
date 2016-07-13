## abhasi, cgates
## 6/2016
## Watermelon: Recreate Legacy pipeline in snakemake

from collections import defaultdict

configfile: "config.yaml"

rule all:
    input:
        expand("01-raw_reads/{sample}_R1.fastq.gz",
                sample=config["samples"]),
        expand("02-fastqc_reads/{sample}_R1_fastqc.html",
                sample=config["samples"]),
#         expand("03-tophat/{sample}/{sample}_accepted_hits.bam",
#                 sample=config["samples"]),
        expand("04-fastqc_align/{sample}_accepted_hits_fastqc.html",
                sample=config["samples"]),
        "05-qc_metrics/Align_summary_all.txt",
#        expand("06-htseq/{sample}_counts.txt",
#                sample=config["samples"]),
        "06-htseq/HTSeq_counts.txt",
#         expand("07-cuffdiff/{comparison}/gene_exp.diff",
#                 comparison=config["comparisons"]),
#         expand("07-cuffdiff/{comparison}/isoform_exp.diff",
#                 comparison=config["comparisons"]),
        expand("07-cuffdiff/{comparison}/{comparison}_gene.foldchange.{fold_change}.txt",
                comparison=config["comparisons"],
                fold_change=config["fold_change"]),
        expand("07-cuffdiff/{comparison}/{comparison}_isoform.foldchange.{fold_change}.txt",
                comparison=config["comparisons"],
                fold_change=config["fold_change"])
        # expand("08-gene_annotation/{comparison}/{comparison}_DEG_ids.txt",
        #         comparison=config["comparisons"]),
        # expand("08-gene_annotation/{comparison}/{comparison}_DEG_names.txt",
        #         comparison=config["comparisons"])
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
        "03-tophat"
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

# samples:
#     Sample_01: Wt
#     Sample_02: Mut1
#     Sample_03: Mut2
#     Sample_04: Mut2

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
#        "07-cuffdiff/{comparison}/read_groups.info"
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
        "07-cuffdiff/{comparison}/{comparison}_gene.foldchange.{fold_change}.txt",
        "07-cuffdiff/{comparison}/{comparison}_isoform.foldchange.{fold_change}.txt"
    params:
        output_dir = "07-cuffdiff/{comparison}",
        fold_change = config["fold_change"],
        comparison = lambda wildcards: wildcards.comparison
    shell:
        " module load rnaseq && "
        "perl scripts/Cuffdiff_out_format_v6.pl "
        " {input.cuffdiff_gene_exp} "
        " {params.comparison}_gene "
        " {params.fold_change} "
        " {params.output_dir} && "
        "perl scripts/Cuffdiff_out_format_v6.pl "
        " {input.cuffdiff_isoform_exp} "
        " {params.comparison}_isoform "
        " {params.fold_change} "
        " {params.output_dir} "


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
