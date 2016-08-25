## abhasi, cgates
## 7/26/2016
## Watermelon: Recreate Legacy pipeline in snakemake

## snakemake --snakefile <snakefile> --configfile <config.yaml> --cores
## snakemake --snakefile Snakefile --configfile tronson_config.yaml  --cores 40 -T -D >workflow_summary.xls

from collections import defaultdict
import csv
import os


rule all:
    input:
        expand("03-fastqc_reads/{sample}_trimmed_R1_fastqc.html",
                sample=config["samples"]),
        expand("05-fastqc_align/{sample}_accepted_hits_fastqc.html",
                sample=config["samples"]),
        "06-qc_metrics/alignment_stats.txt",
        expand("09-diff_expression/{comparison}/{comparison}_gene.foldchange.{fold_change}.txt",
                comparison=config["comparisons"],
                fold_change=config["fold_change"]),
        expand("09-diff_expression/{comparison}/{comparison}_isoform.foldchange.{fold_change}.txt",
                comparison=config["comparisons"],
                fold_change=config["fold_change"]),
        expand("11-cummerbund/{comparison}/Plots/{comparison}_MDSRep.pdf",
                comparison=config["comparisons"]),
        expand("12-deseq2/{comparison}/DESeq2_{comparison}_DE.txt",
                comparison=config["comparisons"]),
        expand("12-deseq2/{comparison}/DESeq2_{comparison}_DE.txt",
                comparison=config["comparisons"])

rule concat_reads:
    input:
        config["input_dir"] + "/{sample}/"
    output:
        "01-raw_reads/{sample}_R1.fastq.gz"
    shell:
        "cat {input}/* > {output}"

rule cutadapt:
    input:
         "01-raw_reads/{sample}_R1.fastq.gz"
    output:
        "02-cutadapt/{sample}_trimmed_R1.fastq.gz"
    params:
        base_quality_5prime_3prime = config["trimming_options"]["base_quality_5prime_3prime"],
        trim_length_5prime = config["trimming_options"]["trim_length_5prime"],
        trim_length_3prime = config["trimming_options"]["trim_length_3prime"]
    log:
        "02-cutadapt/{sample}_trimmed_R1.log"
    shell:
        " module load rnaseq && "
        "cutadapt -q {params.base_quality_5prime_3prime} "
        " -u {params.trim_length_5prime} "
        " -u -{params.trim_length_3prime} "
        " --trim-n -m 20 "
        " -o {output} "
        " {input} "
        " 2>&1 | tee {log}"

rule fastqc_reads:
    input:
        "02-cutadapt/{sample}_trimmed_R1.fastq.gz"
    output:
        "03-fastqc_reads/{sample}_trimmed_R1_fastqc.html"
    log:
        "03-fastqc_reads/{sample}_trimmed_R1_fastqc.log"
    shell:
        " module load rnaseq && "
        "fastqc {input} -o 03-fastqc_reads 2>&1 | tee {log}"

def tophat_options(alignment_options):
    options = ''
    if not isinstance(alignment_options['transcriptome_only'], bool):
        raise ValueError("config alignment_options:transcriptome_only must be boolean")
    if alignment_options['transcriptome_only']:
        options += ' --transcriptome-only '
    return options

rule tophat:
    input:
        "02-cutadapt/{sample}_trimmed_R1.fastq.gz"
    output:
        "04-tophat/{sample}/{sample}_accepted_hits.bam",
        "04-tophat/{sample}/{sample}_align_summary.txt"
    params:
        gtf_file = config["gtf"],
        transcriptome_index= config["transcriptome_index"],
        bowtie2_index = config["bowtie2_index"],
        sample = lambda wildcards: wildcards.sample,
        tophat_options = lambda wildcards: tophat_options(config["alignment_options"])
    shell: " module load rnaseq && "
            "tophat {params.tophat_options} "
            " -p 8 "
            " --b2-very-sensitive "
            " --no-coverage-search "
            " --library-type fr-unstranded "
            " -I 500000 "
#            " -G {params.gtf_file} "   #this option will lead to recreation of the index every time; once a transcriptome index is created, don't give -G
            " --transcriptome-index={params.transcriptome_index} "
            " --no-novel-juncs "
            " -o 04-tophat/{params.sample} "
            " {params.bowtie2_index} "
            " {input} && "
            " mv 04-tophat/{params.sample}/accepted_hits.bam 04-tophat/{params.sample}/{params.sample}_accepted_hits.bam && "
            " mv 04-tophat/{params.sample}/align_summary.txt 04-tophat/{params.sample}/{params.sample}_align_summary.txt "

rule fastqc_align:
    input:
        "04-tophat/{sample}/{sample}_accepted_hits.bam"
    output:
        "05-fastqc_align/{sample}_accepted_hits_fastqc.html"
    shell:
        " module load rnaseq && "
        "fastqc {input} -o 05-fastqc_align"

rule align_qc_metrics:
    input:
        expand("04-tophat/{sample}/{sample}_align_summary.txt", sample=config["samples"])
    output:
        "06-qc_metrics/alignment_stats.txt"
    shell:
        "find 04-tophat -name '*align_summary.txt' | "
        "sort | xargs awk "
        "'BEGIN {{print \"sample\tinput_reads\tmapped_reads\talignment_rate\"}} "
        "/Reads/ {{n=split(FILENAME, fields, /\//); printf \"%s\t\",fields[n-1]}} "
        "/Input/ {{printf \"%s\t\",$3}} "
        "/Mapped/ {{printf \"%s\t\",$3}} "
        "/overall/ {{print $1}}' > {output}"

rule htseq_per_sample:
    input:
        "04-tophat/{sample}/{sample}_accepted_hits.bam"
    output:
        "07-htseq/{sample}_counts.txt",
    params:
        gtf = config["gtf"],
        input_dir = "07-htseq",
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
        "perl /ccmb/BioinfCore/SoftwareDev/projects/Watermelon/scripts/mergeHTSeqCountFiles.pl {params.input_dir}"

# unable to do this in one step; gives error: 'Not all output files of rule htseq_per_sample contain the same wildcards.'
rule htseq_merge:
    input:
        expand("07-htseq/{sample}_counts.txt", sample=config["samples"])
    output:
        "07-htseq/HTSeq_counts.txt"
    params:
        output_dir = "07-htseq",
        input_dir = "07-htseq"
    shell:
       " perl /ccmb/BioinfCore/SoftwareDev/projects/Watermelon/scripts/mergeHTSeqCountFiles.pl {params.input_dir} "

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
        expand("04-tophat/{sample}/{sample}_accepted_hits.bam", sample=config["samples"])
    output:
        "08-cuffdiff/{comparison}/gene_exp.diff",
        "08-cuffdiff/{comparison}/isoform_exp.diff",
        "08-cuffdiff/{comparison}/read_groups.info"
    params:
        output_dir = "08-cuffdiff/{comparison}",
        labels = lambda wildcards : cuffdiff_labels(wildcards.comparison),
        samples = lambda wildcards : cuffdiff_samples(wildcards.comparison,
                                                      config["samples"],
                                                      "04-tophat/{sample_placeholder}/{sample_placeholder}_accepted_hits.bam"),
        fasta = config["fasta"],
        gtf_file = config["gtf"],
    shell:
        " module load rnaseq && "
        " cuffdiff -q "
        " -L {params.labels} "
        " --max-bundle-frags 999999999 "
        " -o {params.output_dir} "
        " -b {params.fasta} "
        " -u -N "
        " --compatible-hits-norm "
        " {params.gtf_file} "
        " {params.samples}" 

rule diff_exp:
    input:
        cuffdiff_gene_exp="08-cuffdiff/{comparison}/gene_exp.diff",
        cuffdiff_isoform_exp="08-cuffdiff/{comparison}/isoform_exp.diff"
    output:
        "09-diff_expression/{comparison}/{comparison}_gene.foldchange.{fold_change}.txt",
        "09-diff_expression/{comparison}/{comparison}_isoform.foldchange.{fold_change}.txt"
    params:
        output_dir = "09-diff_expression/{comparison}",
        fold_change = config["fold_change"],
        comparison = lambda wildcards: wildcards.comparison
    shell:
        " module load rnaseq && "
        " mkdir -p {params.output_dir} &&"
        " perl /ccmb/BioinfCore/SoftwareDev/projects/Watermelon/scripts/Cuffdiff_out_format_v6.pl "
        " {input.cuffdiff_gene_exp} "
        " {params.comparison}_gene "
        " {params.fold_change} "
        " {params.output_dir} && "
        " perl /ccmb/BioinfCore/SoftwareDev/projects/Watermelon/scripts/Cuffdiff_out_format_v6.pl "
        " {input.cuffdiff_isoform_exp} "
        " {params.comparison}_isoform "
        " {params.fold_change} "
        " {params.output_dir} "

rule deseq_build_group_replicates:
    input:
        "08-cuffdiff/{comparison}/read_groups.info"
    output:
        "10-group_replicates/{comparison}/group_replicates.txt"
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
        "10-group_replicates/{comparison}/group_replicates.txt"
    output:
        "11-cummerbund/{comparison}/Plots/{comparison}_MDSRep.pdf"
    params:
        cuff_diff_dir = "08-cuffdiff/{comparison}",
        output_dir = "11-cummerbund/{comparison}/",
        gtf_file = config["gtf"],
        genome = config["genome"],
        logfile = "11-cummerbund/{comparison}/cummerbund.log"
    log:
         "11-cummerbund/{comparison}/cummerbund.log"
    shell:
        " module load rnaseq && "
        " mkdir -p {params.output_dir}/Plots && "
        " Rscript /ccmb/BioinfCore/SoftwareDev/projects/Watermelon/scripts/Run_cummeRbund.R "
        " baseDir={params.output_dir} "
        " cuffDiffDir={params.cuff_diff_dir} "
        " grpRepFile={input} "
        " gtfFile={params.gtf_file} "
        " genome={params.genome} "
        " 2> {log} "

rule deseq2:
    input:
        counts_file = "07-htseq/HTSeq_counts.txt",
        group_replicates = "10-group_replicates/{comparison}/group_replicates.txt"
    output:
        "12-deseq2/{comparison}/DESeq2_{comparison}_DE.txt",
        "12-deseq2/{comparison}/DESeq2_{comparison}_DESig.txt"
    params:
        output_dir = "12-deseq2",
        comparison = lambda wildcards: wildcards.comparison
    log:
         "12-deseq2/{comparison}/DESeq2_{comparison}.log"
    shell:
        " module load rnaseq && "
        " Rscript /ccmb/BioinfCore/SoftwareDev/projects/Watermelon/scripts/Run_DESeq.R "
        " baseDir={params.output_dir}/{params.comparison} "
        " grpRepFile={input.group_replicates} "
        " countsFile={input.counts_file} "
        " 2> {log} "
