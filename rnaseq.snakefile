## abhasi, cgates
## 7/26/2016
## Watermelon 1.0 : Recreate Legacy pipeline in snakemake

## snakemake --snakefile <snakefile> --configfile <config.yaml> --cores 40 -T
## snakemake --snakefile Snakefile --configfile tronson_config.yaml  --cores 40 -T -D >workflow_summary.xls
from __future__ import print_function, absolute_import, division

from collections import defaultdict
from itertools import combinations
import csv
import os

import scripts.rnaseq_snakefile_helper as rnaseq_snakefile_helper

WATERMELON_SCRIPTS_DIR = os.environ.get('WATERMELON_SCRIPTS_DIR', 'scripts')
USER_EMAIL = os.getlogin() + '@umich.edu'
PROJECT_DIR = os.path.basename(os.getcwd())

COMPARISON_INFIX = '_v_'

rnaseq_snakefile_helper.checksum_reset_all("config_checksums", config)
rnaseq_snakefile_helper.init_references(config["references"])

rule all:
    input:
        expand("03-fastqc_reads/{sample}_trimmed_R1_fastqc.html",
                sample=config["samples"]),
        expand("05-fastqc_align/{sample}_accepted_hits_fastqc.html",
                sample=config["samples"]),
        "06-qc_metrics/alignment_stats.txt",
        expand("08-cuffdiff/{multi_group_comparison}/gene_exp.diff",
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        expand("10-flag_diff_expression/{multi_group_comparison}/{multi_group_comparison}_gene.flagged.txt",
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        expand("10-flag_diff_expression/{multi_group_comparison}/{multi_group_comparison}_isoform.flagged.txt",
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        expand("11-annotated_flag_diff_expression/{multi_group_comparison}/{multi_group_comparison}_gene.flagged.annot.txt",
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"]),
                fold_change=config["fold_change"]),
        expand("11-annotated_flag_diff_expression/{multi_group_comparison}/{multi_group_comparison}_isoform.flagged.annot.txt",
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"]),
                fold_change=config["fold_change"]),
        expand("13-cummerbund/{multi_group_comparison}/Plots/{multi_group_comparison}_boxplot.pdf",
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        "07-htseq/HTSeq_counts.txt",
        expand("14-split_diff_expression/{user_specified_comparison}_gene.txt", user_specified_comparison = config["comparisons"]),
        expand("14-split_diff_expression/{user_specified_comparison}_isoform.txt", user_specified_comparison = config["comparisons"])

rule concat_reads:
    input:
        config["input_dir"] + "/{sample}/"
    output:
        "01-raw_reads/{sample}_R1.fastq.gz"
    shell:
        "cat {input}/*.fastq.gz > {output}"

rule cutadapt:
    input:
        trimming_options_checksum = "config_checksums/trimming_options.watermelon.md5",
        raw_fastq = "01-raw_reads/{sample}_R1.fastq.gz"
    output:
        "02-cutadapt/{sample}_trimmed_R1.fastq.gz"
    params:
        base_quality_5prime = config["trimming_options"]["base_quality_5prime"],
        base_quality_3prime = config["trimming_options"]["base_quality_3prime"],
        trim_length_5prime = config["trimming_options"]["trim_length_5prime"],
        trim_length_3prime = config["trimming_options"]["trim_length_3prime"],
        output_dir = "02-cutadapt",
        output_file = "{sample}_trimmed_R1.fastq.gz",
        trimming_options = rnaseq_snakefile_helper.cutadapt_options(config["trimming_options"])
    log:
        "02-cutadapt/{sample}_cutadapt.log"
    shell:
        "module load rnaseq && "
        "if [[ {params.trimming_options} > 0 ]]; then "
        "(cutadapt -q {params.base_quality_5prime},{params.base_quality_3prime} "
        " -u {params.trim_length_5prime} "
        " -u -{params.trim_length_3prime} "
        " --trim-n -m 20 "
        " -o {output}.tmp.gz "
        " {input.raw_fastq} && "
        " mv {output}.tmp.gz {output}) 2>&1 | tee {log}; "
        " else "
        "cd {params.output_dir}; "
        " ln -sf ../{input.raw_fastq} {params.output_file}; "
        " cd ..;"
        " echo \"No trimming done\" 2>&1 |tee {log}; "
        " fi"

rule fastqc_trimmed_reads:
    input:
        "02-cutadapt/{sample}_trimmed_R1.fastq.gz"
    output:
        "03-fastqc_reads/{sample}_trimmed_R1_fastqc.html"
    log:
        "03-fastqc_reads/{sample}_fastqc_trimmed_reads.log"
    shell:
        " module load rnaseq && "
        "fastqc {input} -o 03-fastqc_reads 2>&1 | tee {log}"

rule create_transcriptome_index:
    input:
        alignment_options_checksum = "config_checksums/alignment_options.watermelon.md5",
        reference_checksum = "config_checksums/references.watermelon.md5",
        gtf = "references/gtf",
        bowtie2_index_dir = "references/bowtie2_index"
    output:
        "04-tophat/transcriptome_index/transcriptome.fa"
    params:
        transcriptome_dir = "transcriptome_index",
        temp_dir =  "04-tophat/.tmp",
        output_dir = "04-tophat",
        strand = rnaseq_snakefile_helper.check_strand_option("tuxedo", config["alignment_options"]["library_type"]) 
    log: 
        "04-tophat/create_transcriptome_index.log"
    shell:
        "mkdir -p {params.temp_dir} && "
        " rm -rf {params.temp_dir}/* && "
        " module load rnaseq && "
        " tophat -G {input.gtf} "
        " --library-type {params.strand} "
        " --transcriptome-index={params.temp_dir}/transcriptome_index/transcriptome "
        " {input.bowtie2_index_dir}/genome  2>&1 | tee {log} && "
        " rm -rf {params.output_dir}/{params.transcriptome_dir} && "
        " mv {params.temp_dir}/{params.transcriptome_dir} {params.output_dir} && "
        " mv tophat_out {params.output_dir}/transcriptome_index/ && "
        " touch {params.output_dir}/transcriptome_index/* "

rule tophat:
    input:
        alignment_options_checksum = "config_checksums/alignment_options.watermelon.md5",
        reference_checksum = "config_checksums/references.watermelon.md5",
        transcriptome_fasta = "04-tophat/transcriptome_index/transcriptome.fa",
        bowtie2_index_dir = "references/bowtie2_index",
        fastq = "02-cutadapt/{sample}_trimmed_R1.fastq.gz"
    output:
        "04-tophat/{sample}/{sample}_accepted_hits.bam",
        "04-tophat/{sample}/{sample}_align_summary.txt"
    params:
        transcriptome_index = "04-tophat/transcriptome_index/transcriptome",
        sample = lambda wildcards: wildcards.sample,
        tophat_options = lambda wildcards: rnaseq_snakefile_helper.tophat_options(config["alignment_options"]),
        strand = rnaseq_snakefile_helper.check_strand_option("tuxedo", config["alignment_options"]["library_type"])
    log:
        "04-tophat/{sample}/{sample}_tophat.log"
    threads: 8
    shell: 
            " module load rnaseq && "
            " tophat -p {threads} "
            " --b2-very-sensitive "
            " --no-coverage-search "
            " --library-type {params.strand} "
            " -I 500000 "
            " --transcriptome-index={params.transcriptome_index} "
            " {params.tophat_options} "
            " -o 04-tophat/{params.sample} "
            " {input.bowtie2_index_dir}/genome "
            " {input.fastq} "
            " 2>&1 | tee {log} && "
            " mv 04-tophat/{params.sample}/accepted_hits.bam 04-tophat/{params.sample}/{params.sample}_accepted_hits.bam && "
            " mv 04-tophat/{params.sample}/align_summary.txt 04-tophat/{params.sample}/{params.sample}_align_summary.txt "

rule fastqc_tophat_align:
    input:
        "04-tophat/{sample}/{sample}_accepted_hits.bam"
    output:
        "05-fastqc_align/{sample}_accepted_hits_fastqc.html"
    log:
        "05-fastqc_align/{sample}_fastqc_tophat_align.log"
    shell:
        " module load rnaseq && "
        "fastqc {input} -o 05-fastqc_align 2>&1 | tee {log} "

rule align_qc_metrics:
    input:
        sample_checksum = "config_checksums/samples.watermelon.md5",
        align_summary_files = expand("04-tophat/{sample}/{sample}_align_summary.txt", sample=config["samples"])
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
        reference_checksum = "config_checksums/references.watermelon.md5",
        bams = "04-tophat/{sample}/{sample}_accepted_hits.bam",
        gtf = "references/gtf"
    output:
        "07-htseq/{sample}_counts.txt"
    params:
        input_dir = "07-htseq",
        strand = rnaseq_snakefile_helper.check_strand_option("htseq", config["alignment_options"]["library_type"])
    log:
        "07-htseq/{sample}_htseq_per_sample.log"
    shell:
        " module load rnaseq &&"
        " python -m HTSeq.scripts.count "
        " -f bam "
        " -s {params.strand} "
        " -m "
        " intersection-nonempty "
        " -q {input.bams} "
        " {input.gtf} "
        " > {output} "
        " 2> {log} "

rule htseq_merge:
    input:
        sample_checksum = "config_checksums/samples.watermelon.md5",
        sample_count_files = expand("07-htseq/{sample}_counts.txt", sample=config["samples"])
    output:
        "07-htseq/HTSeq_counts.txt"
    params:
        output_dir = "07-htseq",
        input_dir = "07-htseq"
    shell:
       " perl {WATERMELON_SCRIPTS_DIR}/mergeHTSeqCountFiles.pl {params.input_dir} "

rule cuffdiff:
    input:
        sample_checksum = "config_checksums/samples.watermelon.md5",
        comparison_checksum = "config_checksums/comparisons.watermelon.md5",
        reference_checksum = "config_checksums/references.watermelon.md5",
        fasta_file = "references/bowtie2_index/genome.fa",
        gtf_file = "references/gtf",
        bam_files = expand("04-tophat/{sample}/{sample}_accepted_hits.bam", sample=config["samples"])
    output:
        "08-cuffdiff/{multi_group_comparison}/gene_exp.diff",
        "08-cuffdiff/{multi_group_comparison}/isoform_exp.diff",
        "08-cuffdiff/{multi_group_comparison}/read_groups.info"
    params:
        output_dir = "08-cuffdiff/{multi_group_comparison}",
        labels = lambda wildcards : rnaseq_snakefile_helper.cuffdiff_labels(COMPARISON_INFIX, wildcards.multi_group_comparison),
        samples = lambda wildcards : rnaseq_snakefile_helper.cuffdiff_samples(COMPARISON_INFIX,
                                                                              wildcards.multi_group_comparison,
                                                                              config["samples"],
                                                                              "04-tophat/{sample_placeholder}/{sample_placeholder}_accepted_hits.bam"),
        strand = rnaseq_snakefile_helper.check_strand_option("tuxedo", config["alignment_options"]["library_type"])
    threads: 8
    log:
        "08-cuffdiff/{multi_group_comparison}/{multi_group_comparison}_cuffdiff.log"
    shell:
        " module load rnaseq && "
        " cuffdiff -q "
        " -p {threads} "
        " -L {params.labels} "
        " --max-bundle-frags 999999999 "
        " --library-type {params.strand} "
        " -o {params.output_dir} "
        " -b {input.fasta_file} "
        " -u -N "
        " --compatible-hits-norm "
        " {input.gtf_file} "
        " {params.samples} "
        " 2>&1 | tee {log} "

rule flip_diffex:
    input:
        gene_cuffdiff = "08-cuffdiff/{multi_group_comparison}/gene_exp.diff",
        isoform_cuffdiff = "08-cuffdiff/{multi_group_comparison}/isoform_exp.diff"
    output:
        gene_flip = "09-flip_diffex/{multi_group_comparison}/gene_exp.flip.diff",
        isoform_flip = "09-flip_diffex/{multi_group_comparison}/isoform_exp.flip.diff"
    params:
        comparisons = ",".join(config["comparisons"].values())
    shell:
        " module purge && module load python/3.4.3 && "
        "python {WATERMELON_SCRIPTS_DIR}/flip_diffex.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " {input.gene_cuffdiff} "
        " {output.gene_flip} "
        " {params.comparisons} && "
        
        "python {WATERMELON_SCRIPTS_DIR}/flip_diffex.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " {input.isoform_cuffdiff} "
        " {output.isoform_flip} "
        " {params.comparisons} "

rule flag_diffex:
    input:
        fold_change_checksum = "config_checksums/fold_change.watermelon.md5",
        cuffdiff_gene_exp = "09-flip_diffex/{comparison}/gene_exp.flip.diff",
        cuffdiff_isoform_exp = "09-flip_diffex/{comparison}/isoform_exp.flip.diff"
    output:
        gene_flagged = "10-flag_diff_expression/{comparison}/{comparison}_gene.flagged.txt",
        isoform_flagged = "10-flag_diff_expression/{comparison}/{comparison}_isoform.flagged.txt",
    params:
        fold_change = config["fold_change"]
    log:
        "10-flag_diff_expression/{comparison}/{comparison}_flag_diffex.log"
    shell: 
        " module purge && "
        " module load python/3.4.3 && "
        " python {WATERMELON_SCRIPTS_DIR}/flag_diffex.py "
        " -f {params.fold_change} "
        " {input.cuffdiff_gene_exp} "
        " {output.gene_flagged} "
        " 2>&1 | tee {log} && "
        
        " python {WATERMELON_SCRIPTS_DIR}/flag_diffex.py "
        " -f {params.fold_change} "
        " {input.cuffdiff_isoform_exp} "
        " {output.isoform_flagged} "
        " 2>&1 | tee >>{log} "

rule annotate_flag_diffex:
    input:
        genome_checksum = "config_checksums/genome.watermelon.md5",
        reference_checksum = "config_checksums/references.watermelon.md5",
        gene_diff_exp = "10-flag_diff_expression/{comparison}/{comparison}_gene.flagged.txt",
        isoform_diff_exp = "10-flag_diff_expression/{comparison}/{comparison}_isoform.flagged.txt",
        entrez_gene_info = "references/entrez_gene_info"
    output:
        gene_annot = "11-annotated_flag_diff_expression/{comparison}/{comparison}_gene.flagged.annot.txt",
        isoform_annot = "11-annotated_flag_diff_expression/{comparison}/{comparison}_isoform.flagged.annot.txt"
    params:
        output_dir = "11-annotated_flag_diff_expression/{comparison}",
        genome = config["genome"]
    log:
        "11-annotated_flag_diff_expression/{comparison}/{comparison}_annotate_flag_diffex.log"
    shell:
        "python {WATERMELON_SCRIPTS_DIR}/annotate_entrez_gene_info.py "
        " -i {input.entrez_gene_info} "
        " -e {input.gene_diff_exp} "
        " -g {params.genome} "
        " -o {params.output_dir} "
        " 2>&1 | tee {log} && "
        
        " python {WATERMELON_SCRIPTS_DIR}/annotate_entrez_gene_info.py "
        " -i {input.entrez_gene_info} "
        " -e {input.isoform_diff_exp} "
        " -g {params.genome} "
        " -o {params.output_dir} "
        " 2>&1 | tee >>{log} "

rule build_group_replicates:
    input:
        "08-cuffdiff/{comparison}/read_groups.info"
    output:
        "12-group_replicates/{comparison}/group_replicates.txt"
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
        genome_checksum = "config_checksums/genome.watermelon.md5",
        reference_checksum = "config_checksums/references.watermelon.md5",
        group_replicates = "12-group_replicates/{comparison}/group_replicates.txt",
        gtf_file = "references/gtf"
    output:
        "13-cummerbund/{comparison}/Plots/{comparison}_boxplot.pdf" #should we list out all outputs here?
    params:
        cuff_diff_dir = "08-cuffdiff/{comparison}",
        output_dir = "13-cummerbund/{comparison}/",
        genome = config["genome"],
        logfile = "13-cummerbund/{comparison}/cummerbund.log"
    log:
         "13-cummerbund/{comparison}/{comparison}_cummerbund.log"
    shell:
        "module load rnaseq && "
        "mkdir -p {params.output_dir}/Plots && "
        "Rscript {WATERMELON_SCRIPTS_DIR}/Run_cummeRbund.R "
        " baseDir={params.output_dir} "
        " cuffDiffDir={params.cuff_diff_dir} "
        " grpRepFile={input.group_replicates} "
        " gtfFile={input.gtf_file} "
        " genome={params.genome} "
        " 2>&1 | tee {log} "

rule split_diffex:
    input:
        gene = expand("11-annotated_flag_diff_expression/{multi_group_comparison}/{multi_group_comparison}_gene.flagged.annot.txt",
                            multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        isoform =  expand("11-annotated_flag_diff_expression/{multi_group_comparison}/{multi_group_comparison}_isoform.flagged.annot.txt",
                        multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"]))
    output:
        expand("14-split_diff_expression/{user_specified_comparisons}_gene.txt", user_specified_comparisons=config["comparisons"]),
        expand("14-split_diff_expression/{user_specified_comparisons}_isoform.txt", user_specified_comparisons=config["comparisons"]),
    params:
        output_dir = "14-split_diff_expression",
        user_specified_comparison_list = ",".join(config["comparisons"].values())
    shell:
        "module purge && module load python/3.4.3 && "
        "python {WATERMELON_SCRIPTS_DIR}/split_diffex.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " -o _gene.txt "
        " {input.gene} "
        " {params.output_dir} "
        " {params.user_specified_comparison_list} && "
        
        "module purge && module load python/3.4.3 && "
        "python {WATERMELON_SCRIPTS_DIR}/split_diffex.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " -o _isoform.txt "
        " {input.isoform} "
        " {params.output_dir} "
        " {params.user_specified_comparison_list} "

onsuccess:
    print("Workflow finished, no errors")
    shell("echo 'Run complete for {PROJECT_DIR}. No errors'  | mail -s 'Project {PROJECT_DIR}: Run completed for Watermelon RNA-Seq pipeline' {USER_EMAIL}")
onerror:
    print("Workflow incomplete: Error occurred. See watermelon.log for details")
    shell("echo 'Watermelon encountered an unexpected error for {PROJECT_DIR}. Review config and contact bfxcore support' | mail -s 'Error in Watermelon RNA-Seq pipeline for {PROJECT_DIR}' {USER_EMAIL}")

