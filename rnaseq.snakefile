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
import subprocess
import yaml

import scripts.rnaseq_snakefile_helper as rnaseq_snakefile_helper

WATERMELON_SCRIPTS_DIR = os.environ.get('WATERMELON_SCRIPTS_DIR', 'scripts')

COMPARISON_INFIX = '_v_'
INPUT_DIR = config.get("input_dir", "inputs")
ALIGNMENT_DIR = config.get("alignment_output_dir", "alignment_results")
DIFFEX_DIR = config.get("diffex_output_dir", "diffex_results")

rnaseq_snakefile_helper.checksum_reset_all("config_checksums", config)
rnaseq_snakefile_helper.init_references(config["references"])

rule all:
    input:
        expand("{alignment_dir}/03-fastqc_reads/{sample}_trimmed_R1_fastqc.html",
                alignment_dir=ALIGNMENT_DIR,
                sample=config["samples"]),
        expand("{alignment_dir}/04-tophat/{sample}/{sample}_accepted_hits.bam",
                alignment_dir=ALIGNMENT_DIR,
                sample=config["samples"]),
        expand("{alignment_dir}/05-fastqc_align/{sample}_accepted_hits_fastqc.html",
                alignment_dir=ALIGNMENT_DIR,
                sample=config["samples"]),
        ALIGNMENT_DIR + "/06-qc_metrics/alignment_stats.txt",
        DIFFEX_DIR + "/07-htseq/HTSeq_counts.txt",
        expand("{diffex_dir}/08-cuffdiff/{multi_group_comparison}/gene_exp.diff",
                diffex_dir=DIFFEX_DIR, 
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        expand("{diffex_dir}/10-diffex_flag/{multi_group_comparison}/{multi_group_comparison}_gene.flagged.txt",
                diffex_dir=DIFFEX_DIR, 
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        expand("{diffex_dir}/10-diffex_flag/{multi_group_comparison}/{multi_group_comparison}_isoform.flagged.txt",
                diffex_dir=DIFFEX_DIR,
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        expand("{diffex_dir}/11-annotate_diffex_flag/{multi_group_comparison}/{multi_group_comparison}_gene.flagged.annot.txt",
                diffex_dir=DIFFEX_DIR,
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"]),
                fold_change=config["fold_change"]),
        expand("{diffex_dir}/11-annotate_diffex_flag/{multi_group_comparison}/{multi_group_comparison}_isoform.flagged.annot.txt",
                diffex_dir=DIFFEX_DIR,
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"]),
                fold_change=config["fold_change"]),
        expand("{diffex_dir}/13-cummerbund/{multi_group_comparison}/Plots/{multi_group_comparison}_boxplot.pdf",
                diffex_dir=DIFFEX_DIR,
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        expand("{diffex_dir}/13-cummerbund/{multi_group_comparison}/{multi_group_comparison}_repRawCounts.txt",
                diffex_dir=DIFFEX_DIR,
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        DIFFEX_DIR + "/Deliverables/qc/raw_reads_fastqc",
        DIFFEX_DIR + "/Deliverables/qc/aligned_reads_fastqc",
        DIFFEX_DIR + "/Deliverables/qc/alignment_stats.txt",
        expand("{diffex_dir}/Deliverables/diffex/cuffdiff_results/{user_specified_comparisons}.xlsx", 
                diffex_dir=DIFFEX_DIR,
                user_specified_comparisons=config["comparisons"]),
        expand("{diffex_dir}/Deliverables/diffex/{multi_group_comparison}_repRawCounts.txt", 
                diffex_dir=DIFFEX_DIR,
                multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        DIFFEX_DIR + "/Deliverables/diffex/cummeRbund_plots"


rule concat_reads:
    input:
        INPUT_DIR + "/{sample}/"
    output:
        ALIGNMENT_DIR + "/01-raw_reads/{sample}_R1.fastq.gz"
    shell:
        "cat {input}/*.fastq.gz > {output}"

rule cutadapt:
    input:
        trimming_options_checksum = "config_checksums/trimming_options.watermelon.md5",
        raw_fastq = ALIGNMENT_DIR + "/01-raw_reads/{sample}_R1.fastq.gz"
    output:
        ALIGNMENT_DIR + "/02-cutadapt/{sample}_trimmed_R1.fastq.gz"
    params:
        base_quality_5prime = config["trimming_options"]["base_quality_5prime"],
        base_quality_3prime = config["trimming_options"]["base_quality_3prime"],
        trim_length_5prime = config["trimming_options"]["trim_length_5prime"],
        trim_length_3prime = config["trimming_options"]["trim_length_3prime"],
        output_dir = ALIGNMENT_DIR + "/02-cutadapt",
        output_file = "{sample}_trimmed_R1.fastq.gz",
        trimming_options = rnaseq_snakefile_helper.cutadapt_options(config["trimming_options"])
    log:
        ALIGNMENT_DIR + "/02-cutadapt/{sample}_cutadapt.log"
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
        " ln -sf ../../{input.raw_fastq} {params.output_dir}/{params.output_file}; "
        " echo \"No trimming done\" 2>&1 |tee {log}; "
        " fi"

rule fastqc_trimmed_reads:
    input:
        ALIGNMENT_DIR + "/02-cutadapt/{sample}_trimmed_R1.fastq.gz"
    output:         
        touch(ALIGNMENT_DIR + "/03-fastqc_reads/reads_fastq.done"),
        ALIGNMENT_DIR + "/03-fastqc_reads/{sample}_trimmed_R1_fastqc.html"
    log:
        ALIGNMENT_DIR + "/03-fastqc_reads/{sample}_fastqc_trimmed_reads.log"
    params:
        fastqc_dir = ALIGNMENT_DIR + "/03-fastqc_reads"
    shell:
        " module load rnaseq && "
        "fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log}"

rule create_transcriptome_index:
    input:
        alignment_options_checksum =  "config_checksums/alignment_options.watermelon.md5",
        reference_checksum =  "config_checksums/references.watermelon.md5",
        gtf = "references/gtf",
        bowtie2_index_dir = "references/bowtie2_index"
    output:
        ALIGNMENT_DIR + "/04-tophat/transcriptome_index/transcriptome.fa"
    params:
        transcriptome_dir = "transcriptome_index",
        temp_dir =  ALIGNMENT_DIR + "/04-tophat/.tmp",
        output_dir = ALIGNMENT_DIR + "/04-tophat",
        strand = rnaseq_snakefile_helper.check_strand_option("tuxedo", config["alignment_options"]["library_type"]) 
    log: 
        ALIGNMENT_DIR + "/04-tophat/create_transcriptome_index.log"
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
        " mv tophat_out {params.output_dir}/{params.transcriptome_dir}/ && "
        " touch {params.output_dir}/{params.transcriptome_dir}/* "

rule tophat:
    input:
        alignment_options_checksum = "config_checksums/alignment_options.watermelon.md5",
        reference_checksum = "config_checksums/references.watermelon.md5",
        transcriptome_fasta = ALIGNMENT_DIR + "/04-tophat/transcriptome_index/transcriptome.fa",
        bowtie2_index_dir = "references/bowtie2_index",
        fastq = ALIGNMENT_DIR + "/02-cutadapt/{sample}_trimmed_R1.fastq.gz"
    output:
        ALIGNMENT_DIR + "/04-tophat/{sample}/{sample}_accepted_hits.bam",
        ALIGNMENT_DIR + "/04-tophat/{sample}/{sample}_align_summary.txt"
    params:
        transcriptome_index = ALIGNMENT_DIR + "/04-tophat/transcriptome_index/transcriptome",
        tophat_dir = ALIGNMENT_DIR + "/04-tophat",
        sample = lambda wildcards: wildcards.sample,
        tophat_options = lambda wildcards: rnaseq_snakefile_helper.tophat_options(config["alignment_options"]),
        strand = rnaseq_snakefile_helper.check_strand_option("tuxedo", config["alignment_options"]["library_type"])
    log:
        ALIGNMENT_DIR + "/04-tophat/{sample}/{sample}_tophat.log"
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
            " -o {params.tophat_dir}/{params.sample} "
            " {input.bowtie2_index_dir}/genome "
            " {input.fastq} "
            " 2>&1 | tee {log} && "
            " mv {params.tophat_dir}/{params.sample}/accepted_hits.bam {params.tophat_dir}/{params.sample}/{params.sample}_accepted_hits.bam && "
            " mv {params.tophat_dir}/{params.sample}/align_summary.txt {params.tophat_dir}/{params.sample}/{params.sample}_align_summary.txt "

rule fastqc_tophat_align:
    input:
        ALIGNMENT_DIR + "/04-tophat/{sample}/{sample}_accepted_hits.bam"
    output:
        touch(ALIGNMENT_DIR + "/05-fastqc_align/align_fastq.done"),
        ALIGNMENT_DIR + "/05-fastqc_align/{sample}_accepted_hits_fastqc.html"
    params:
        fastqc_dir =  ALIGNMENT_DIR + "/05-fastqc_align"
    log:
        ALIGNMENT_DIR + "/05-fastqc_align/{sample}_fastqc_tophat_align.log"
    shell:
        " module load rnaseq && "
        "fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log} "

rule align_qc_metrics:
    input:
        sample_checksum = "config_checksums/samples.watermelon.md5",
        align_summary_files = expand("{alignment_dir}/04-tophat/{sample}/{sample}_align_summary.txt", 
                                        alignment_dir= ALIGNMENT_DIR, sample=config["samples"])
    output:
        ALIGNMENT_DIR + "/06-qc_metrics/alignment_stats.txt"
    params:
        tophat_dir = ALIGNMENT_DIR + "/04-tophat"
    shell:
        "find {params.tophat_dir} -name '*align_summary.txt' | "
        "sort | xargs awk "
        "'BEGIN {{print \"sample\tinput_reads\tmapped_reads\talignment_rate\"}} "
        "/Reads/ {{n=split(FILENAME, fields, /\//); printf \"%s\t\",fields[n-1]}} "
        "/Input/ {{printf \"%s\t\",$3}} "
        "/Mapped/ {{printf \"%s\t\",$3}} "
        "/overall/ {{print $1}}' > {output}"

rule htseq_per_sample:
    input:
        reference_checksum = "config_checksums/references.watermelon.md5",
        bams = ALIGNMENT_DIR + "/04-tophat/{sample}/{sample}_accepted_hits.bam",
        gtf = "references/gtf"
    output:
        DIFFEX_DIR + "/07-htseq/{sample}_counts.txt"
    params:
        strand = rnaseq_snakefile_helper.check_strand_option("htseq", config["alignment_options"]["library_type"])
    log:
        DIFFEX_DIR + "/07-htseq/{sample}_htseq_per_sample.log"
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
        sample_count_files = expand("{diffex_dir}/07-htseq/{sample}_counts.txt",
                                    diffex_dir=DIFFEX_DIR, sample=config["samples"])
    output:
        DIFFEX_DIR + "/07-htseq/HTSeq_counts.txt"
    params:
        output_dir = DIFFEX_DIR + "/07-htseq",
        input_dir = DIFFEX_DIR + "/07-htseq"
    shell:
       " perl {WATERMELON_SCRIPTS_DIR}/mergeHTSeqCountFiles.pl {params.input_dir} "

rule cuffdiff:
    input:
        sample_checksum = "config_checksums/samples.watermelon.md5",
        comparison_checksum = "config_checksums/comparisons.watermelon.md5",
        reference_checksum = "config_checksums/references.watermelon.md5",
        fasta_file = "references/bowtie2_index/genome.fa",
        gtf_file = "references/gtf",
        bam_files = expand("{alignment_dir}/04-tophat/{sample}/{sample}_accepted_hits.bam",
                            alignment_dir= ALIGNMENT_DIR, sample=config["samples"])
    output:
        DIFFEX_DIR + "/08-cuffdiff/{multi_group_comparison}/gene_exp.diff",
        DIFFEX_DIR + "/08-cuffdiff/{multi_group_comparison}/isoform_exp.diff",
        DIFFEX_DIR + "/08-cuffdiff/{multi_group_comparison}/read_groups.info"
    params:
        output_dir = DIFFEX_DIR + "/08-cuffdiff/{multi_group_comparison}",
        labels = lambda wildcards : rnaseq_snakefile_helper.cuffdiff_labels(COMPARISON_INFIX, wildcards.multi_group_comparison),
        samples = lambda wildcards : rnaseq_snakefile_helper.cuffdiff_samples(COMPARISON_INFIX,
                                                                              wildcards.multi_group_comparison,
                                                                              config["samples"],
                                                                              ALIGNMENT_DIR + "/04-tophat/{sample_placeholder}/{sample_placeholder}_accepted_hits.bam"),
        strand = rnaseq_snakefile_helper.check_strand_option("tuxedo", config["alignment_options"]["library_type"])
    threads: 8
    log:
        DIFFEX_DIR + "/08-cuffdiff/{multi_group_comparison}/{multi_group_comparison}_cuffdiff.log"
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

rule diffex_flip:
    input:
        gene_cuffdiff = DIFFEX_DIR + "/08-cuffdiff/{multi_group_comparison}/gene_exp.diff",
        isoform_cuffdiff = DIFFEX_DIR + "/08-cuffdiff/{multi_group_comparison}/isoform_exp.diff"
    output:
        gene_flip = DIFFEX_DIR + "/09-diffex_flip/{multi_group_comparison}/gene_exp.flip.diff",
        isoform_flip = DIFFEX_DIR + "/09-diffex_flip/{multi_group_comparison}/isoform_exp.flip.diff"
    params:
        comparisons = ",".join(config["comparisons"].values())
    shell:
        " module purge && module load python/3.4.3 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_flip.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " {input.gene_cuffdiff} "
        " {output.gene_flip} "
        " {params.comparisons} && "
        
        "python {WATERMELON_SCRIPTS_DIR}/diffex_flip.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " {input.isoform_cuffdiff} "
        " {output.isoform_flip} "
        " {params.comparisons} "

rule diffex_flag:
    input:
        fold_change_checksum = "config_checksums/fold_change.watermelon.md5",
        cuffdiff_gene_exp = DIFFEX_DIR + "/09-diffex_flip/{comparison}/gene_exp.flip.diff",
        cuffdiff_isoform_exp = DIFFEX_DIR + "/09-diffex_flip/{comparison}/isoform_exp.flip.diff"
    output:
        gene_flagged = DIFFEX_DIR + "/10-diffex_flag/{comparison}/{comparison}_gene.flagged.txt",
        isoform_flagged = DIFFEX_DIR + "/10-diffex_flag/{comparison}/{comparison}_isoform.flagged.txt",
    params:
        fold_change = config["fold_change"]
    log:
        DIFFEX_DIR + "/10-diffex_flag/{comparison}/{comparison}_diffex_flag.log"
    shell: 
        " module purge && "
        " module load python/3.4.3 && "
        " python {WATERMELON_SCRIPTS_DIR}/diffex_flag.py "
        " -f {params.fold_change} "
        " {input.cuffdiff_gene_exp} "
        " {output.gene_flagged} "
        " 2>&1 | tee {log} && "
        
        " python {WATERMELON_SCRIPTS_DIR}/diffex_flag.py "
        " -f {params.fold_change} "
        " {input.cuffdiff_isoform_exp} "
        " {output.isoform_flagged} "
        " 2>&1 | tee >>{log} "

rule annotate_diffex_flag:
    input:
        genome_checksum = "config_checksums/genome.watermelon.md5",
        reference_checksum = "config_checksums/references.watermelon.md5",
        gene_diff_exp = DIFFEX_DIR + "/10-diffex_flag/{comparison}/{comparison}_gene.flagged.txt",
        isoform_diff_exp = DIFFEX_DIR + "/10-diffex_flag/{comparison}/{comparison}_isoform.flagged.txt",
        entrez_gene_info = "references/entrez_gene_info"
    output:
        gene_annot = DIFFEX_DIR + "/11-annotate_diffex_flag/{comparison}/{comparison}_gene.flagged.annot.txt",
        isoform_annot = DIFFEX_DIR + "/11-annotate_diffex_flag/{comparison}/{comparison}_isoform.flagged.annot.txt"
    params:
        output_dir = DIFFEX_DIR + "/11-annotate_diffex_flag/{comparison}",
        genome = config["genome"]
    log:
        DIFFEX_DIR + "/11-annotate_diffex_flag/{comparison}/{comparison}_annotate_diffex_flag.log"
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
        DIFFEX_DIR + "/08-cuffdiff/{comparison}/read_groups.info"
    output:
        DIFFEX_DIR + "/12-group_replicates/{comparison}/group_replicates.txt"
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
        group_replicates = DIFFEX_DIR + "/12-group_replicates/{comparison}/group_replicates.txt",
        gtf_file = "references/gtf"
    output:
        DIFFEX_DIR + "/13-cummerbund/{comparison}/Plots",
        DIFFEX_DIR + "/13-cummerbund/{comparison}/Plots/{comparison}_boxplot.pdf",
        DIFFEX_DIR + "/13-cummerbund/{comparison}/{comparison}_repRawCounts.txt"
    params:
        cuff_diff_dir = DIFFEX_DIR + "/08-cuffdiff/{comparison}",
        output_dir = DIFFEX_DIR + "/13-cummerbund/{comparison}/",
        genome = config["genome"],
    log:
         DIFFEX_DIR + "/13-cummerbund/{comparison}/{comparison}_cummerbund.log"
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

rule diffex_split:
    input:
        gene = expand("{diffex_dir}/11-annotate_diffex_flag/{multi_group_comparison}/{multi_group_comparison}_gene.flagged.annot.txt",
                            diffex_dir=DIFFEX_DIR,
                            multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"])),
        isoform =  expand("{diffex_dir}/11-annotate_diffex_flag/{multi_group_comparison}/{multi_group_comparison}_isoform.flagged.annot.txt",
                        diffex_dir=DIFFEX_DIR,
                        multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, config["comparisons"]))
    output:
        expand("{diffex_dir}/14-diffex_split/{user_specified_comparisons}_gene.txt", 
                diffex_dir=DIFFEX_DIR,
                user_specified_comparisons=config["comparisons"]),
        expand("{diffex_dir}/14-diffex_split/{user_specified_comparisons}_isoform.txt", 
                diffex_dir=DIFFEX_DIR,
                user_specified_comparisons=config["comparisons"]),
        touch(DIFFEX_DIR + "/14-diffex_split/last_split"),
        DIFFEX_DIR + "/14-diffex_split/glossary.txt"
    params:
        output_dir = DIFFEX_DIR + "/14-diffex_split",
        user_specified_comparison_list = ",".join(config["comparisons"].values())
    shell:
        "module purge && module load python/3.4.3 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_split.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " -o _gene.txt "
        " {input.gene} "
        " {params.output_dir} "
        " {params.user_specified_comparison_list} && "
        
        "module purge && module load python/3.4.3 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_split.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " -o _isoform.txt "
        " {input.isoform} "
        " {params.output_dir} "
        " {params.user_specified_comparison_list} && "
        "cp {WATERMELON_SCRIPTS_DIR}/glossary.txt {params.output_dir} "

rule build_run_info:
    input: rules.diffex_split.output
    output: 
        DIFFEX_DIR + "/14-diffex_split/run_info.txt",
    run:
        command = 'module load rnaseq; module list -t 2> {}'.format(output[0])
        subprocess.call(command, shell=True)
        with open(output[0], 'a') as run_info_file: 
            print('\n\nConfig\n', file=run_info_file)
            print(yaml.dump(config, default_flow_style=False), file=run_info_file)

rule diffex_excel:
    input:
        gene = DIFFEX_DIR + "/14-diffex_split/{user_specified_comparisons}_gene.txt",
        isoform = DIFFEX_DIR + "/14-diffex_split/{user_specified_comparisons}_isoform.txt",
        glossary = DIFFEX_DIR + "/14-diffex_split/glossary.txt",
        run_info = DIFFEX_DIR + "/14-diffex_split/run_info.txt"
    output: 
        DIFFEX_DIR + "/15-diffex_excel/{user_specified_comparisons}.xlsx",
    shell:
        "module purge && module load python/3.4.3 && "
        " python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py "
        " -g {input.gene}"
        " -i {input.isoform}"
        " --glossary {input.glossary} "
        " --info_filepath {input.run_info} "
        " {output} "

rule watermelon_deliverables:
    input:
        raw_fastqc = expand("{alignment_dir}/03-fastqc_reads/{sample}_trimmed_R1_fastqc.html",
                                alignment_dir=ALIGNMENT_DIR, sample=config["samples"]),
        align_fastqc = expand("{alignment_dir}/05-fastqc_align/{sample}_accepted_hits_fastqc.html",
                                alignment_dir=ALIGNMENT_DIR, sample=config["samples"]),
        alignment_stats = ALIGNMENT_DIR + "/06-qc_metrics/alignment_stats.txt",
        diffex_excel = expand("{diffex_dir}/15-diffex_excel/{user_specified_comparisons}.xlsx", 
                                diffex_dir=DIFFEX_DIR, user_specified_comparisons=config["comparisons"]),
        diffex_raw_counts = expand("{diffex_dir}/13-cummerbund/{multi_group_comparison}/{multi_group_comparison}_repRawCounts.txt",
                                    diffex_dir=DIFFEX_DIR, 
                                    multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, 
                                                                                                        config["comparisons"])),
        plots = expand("{diffex_dir}/13-cummerbund/{multi_group_comparison}/Plots",
                                    diffex_dir=DIFFEX_DIR, 
                                    multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, 
                                                                                                        config["comparisons"])),
    output:
        raw_fastqc = DIFFEX_DIR + "/Deliverables/qc/raw_reads_fastqc",
        align_fastqc = DIFFEX_DIR + "/Deliverables/qc/aligned_reads_fastqc",
        alignment_stats = DIFFEX_DIR + "/Deliverables/qc/alignment_stats.txt",
        diffex_excel = expand("{diffex_dir}/Deliverables/diffex/cuffdiff_results/{user_specified_comparisons}.xlsx", 
                                diffex_dir=DIFFEX_DIR, 
                                user_specified_comparisons=config["comparisons"]),
        diffex_raw_counts = expand("{diffex_dir}/Deliverables/diffex/{multi_group_comparison}_repRawCounts.txt", 
                                    diffex_dir=DIFFEX_DIR, 
                                    multi_group_comparison=rnaseq_snakefile_helper.cuffdiff_conditions(COMPARISON_INFIX, 
                                                                                                        config["comparisons"])),
        plots = DIFFEX_DIR + "/Deliverables/diffex/cummeRbund_plots"
    params:
        align_fastqc_input_dir = ALIGNMENT_DIR + "/05-fastqc_align",
        raw_fastqc_input_dir = ALIGNMENT_DIR + "/03-fastqc_reads",
        diffex_excel_input_dir = DIFFEX_DIR + "/15-diffex_excel",
        diffex_excel_output_dir = DIFFEX_DIR + "/Deliverables/diffex/cuffdiff_results"
    shell:
        "rm -rf {DIFFEX_DIR}/Deliverables && mkdir -p {DIFFEX_DIR}/Deliverables/qc && "
        " mkdir -p {DIFFEX_DIR}/Deliverables/diffex && "
        " ln -s ../../../{params.raw_fastqc_input_dir} {output.raw_fastqc} && "
        " ln -s ../../../{params.align_fastqc_input_dir} {output.align_fastqc}  && "
        " ln -s ../../../{input.alignment_stats} {output.alignment_stats} && "
        " ln -s ../../../{params.diffex_excel_input_dir} {params.diffex_excel_output_dir} && "
        " ln -s ../../../{input.diffex_raw_counts} {output.diffex_raw_counts} && "
        " ln -s ../../../{input.plots}  {output.plots} "
