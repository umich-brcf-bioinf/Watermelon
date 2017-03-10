## abhasi, cgates
## 7/26/2016
## Watermelon 1.0 : Recreate Legacy pipeline in snakemake

from __future__ import print_function, absolute_import, division

from collections import defaultdict, OrderedDict
from itertools import combinations, repeat
from shutil import copyfile
import csv
import os
import subprocess
import yaml

import scripts.rnaseq_snakefile_helper as rnaseq_snakefile_helper
import scripts.deseq2_helper as deseq2_helper

WATERMELON_SCRIPTS_DIR = os.environ.get('WATERMELON_SCRIPTS_DIR', 'scripts')

COMPARISON_INFIX = '_v_'
INPUT_DIR = config.get("input_dir", "inputs")
ALIGNMENT_DIR = config.get("alignment_output_dir", "alignment_results")
DIFFEX_DIR = config.get("diffex_output_dir", "diffex_results")

DELIMITER = '^'
SAMPLES_KEY = 'samples'
PHENOTYPES_KEY = 'phenotypes'
COMPARISONS_KEY ='comparisons'

phenotypeManager = rnaseq_snakefile_helper.PhenotypeManager(config,
                                                            DELIMITER,
                                                            COMPARISON_INFIX)

PHENOTYPE_NAMES = phenotypeManager.phenotypes_comparisons_tuple.phenotypes
COMPARISON_GROUPS = phenotypeManager.phenotypes_comparisons_tuple.comparisons

rnaseq_snakefile_helper.init_references(config["references"])
rnaseq_snakefile_helper.checksum_reset_all("config_checksums",
                                           config=config,
                                           phenotype_comparisons=config['comparisons'],
                                           phenotype_samples=phenotypeManager.phenotype_sample_list)

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
        expand("{diffex_dir}/08-cuffdiff/{phenotype}/gene_exp.diff",
                diffex_dir=DIFFEX_DIR, 
                phenotype=sorted(config[COMPARISONS_KEY].keys())),
        expand("{diffex_dir}/10-diffex_flag/{phenotype}/{phenotype}_gene.flagged.txt",
                diffex_dir=DIFFEX_DIR,
                phenotype=sorted(config[COMPARISONS_KEY].keys())),
        expand("{diffex_dir}/10-diffex_flag/{phenotype}/{phenotype}_isoform.flagged.txt",
                diffex_dir=DIFFEX_DIR,
                phenotype=sorted(config[COMPARISONS_KEY].keys())),                
        expand("{diffex_dir}/11-annotate_diffex_flag/{phenotype}/{phenotype}_gene.flagged.annot.txt",
                diffex_dir=DIFFEX_DIR,
                phenotype=sorted(config[COMPARISONS_KEY].keys())),
        expand("{diffex_dir}/11-annotate_diffex_flag/{phenotype}/{phenotype}_isoform.flagged.annot.txt",
                diffex_dir=DIFFEX_DIR,
                phenotype=sorted(config[COMPARISONS_KEY].keys())),
        expand("{diffex_dir}/13-cummerbund/{pheno}/Plots/{pheno}_boxplot.pdf",
                diffex_dir=DIFFEX_DIR,
                pheno=sorted(config[COMPARISONS_KEY].keys())),
        expand("{diffex_dir}/13-cummerbund/{pheno}/{pheno}_repRawCounts.txt",
                diffex_dir=DIFFEX_DIR,
                pheno=sorted(config[COMPARISONS_KEY].keys())),
        expand("{diffex_dir}/14-diffex_split/{phenotype_name}/{comparison}_gene.txt",
                zip,
                diffex_dir=repeat(DIFFEX_DIR, len(PHENOTYPE_NAMES)),
                phenotype_name=PHENOTYPE_NAMES, 
                comparison=COMPARISON_GROUPS),
        expand("{diffex_dir}/14-diffex_split/{phenotype_name}/{comparison}_isoform.txt",
                zip,
                diffex_dir=repeat(DIFFEX_DIR, len(PHENOTYPE_NAMES)),
                phenotype_name=PHENOTYPE_NAMES, 
                comparison=COMPARISON_GROUPS),
        DIFFEX_DIR + "/14-diffex_split/last_split",
        expand("{diffex_dir}/16-diffex_excel/{phenotype_name}/{comparison}.xlsx",
                zip,
                diffex_dir=repeat(DIFFEX_DIR, len(PHENOTYPE_NAMES)),
                phenotype_name=PHENOTYPE_NAMES, 
                comparison=COMPARISON_GROUPS),
        expand(DIFFEX_DIR + "/Deliverables/alignment_deliverables/raw_reads_fastqc/{sample}_trimmed_R1_fastqc.html",
                sample=config["samples"]),
        expand(DIFFEX_DIR + "/Deliverables/alignment_deliverables/aligned_reads_fastqc/{sample}_accepted_hits_fastqc.html",
                sample=config["samples"]),
         DIFFEX_DIR + "/Deliverables/alignment_deliverables/alignment_stats.txt",
        expand("{diffex_dir}/Deliverables/diffex_deliverables/cuffdiff_results/{phenotype_name}/{comparison}.xlsx", 
                zip,
                diffex_dir=repeat(DIFFEX_DIR, len(PHENOTYPE_NAMES)),
                phenotype_name=PHENOTYPE_NAMES, comparison=COMPARISON_GROUPS),
        expand("{diffex_dir}/Deliverables/diffex_deliverables/cummeRbund_results/{phenotype_name}_repRawCounts.txt", 
                diffex_dir=DIFFEX_DIR,
                phenotype_name=sorted(config[COMPARISONS_KEY].keys())),
        expand("{diffex_dir}/Deliverables/diffex_deliverables/cummeRbund_results/{phenotype_name}_cummeRbund_plots",
                    diffex_dir=DIFFEX_DIR, 
                    phenotype_name=sorted(config[COMPARISONS_KEY].keys())),
        DIFFEX_DIR + "/deseq2/01-htseq/HTSeq_counts.txt",
        DIFFEX_DIR + "/deseq2/02-metadata_contrasts/sample_metadata.txt",
        DIFFEX_DIR + "/deseq2/02-metadata_contrasts/contrasts.txt"

rule concat_reads:
    input:
        INPUT_DIR + "/{sample}/"
    output:
        ALIGNMENT_DIR + "/01-raw_reads/{sample}_R1.fastq.gz"
    shell:
        "cat {input}/*.fastq.gz > {output}"

rule cutadapt:
    input:
        trimming_options_checksum = "config_checksums/config-trimming_options.watermelon.md5",
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
        ALIGNMENT_DIR + "/02-cutadapt/log/{sample}_cutadapt.log"
    shell:
        "module load watermelon_rnaseq && "
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
        ALIGNMENT_DIR + "/03-fastqc_reads/log/{sample}_fastqc_trimmed_reads.log"
    params:
        fastqc_dir = ALIGNMENT_DIR + "/03-fastqc_reads"
    shell:
        " module load watermelon_rnaseq && "
        "fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log}"

rule create_transcriptome_index:
    input:
        alignment_options_checksum =  "config_checksums/config-alignment_options.watermelon.md5",
        reference_checksum =  "config_checksums/config-references.watermelon.md5",
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
        ALIGNMENT_DIR + "/04-tophat/log/create_transcriptome_index.log"
    shell:
        "mkdir -p {params.temp_dir} && "
        " rm -rf {params.temp_dir}/* && "
        " module load watermelon_rnaseq && "
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
        alignment_options_checksum = "config_checksums/config-alignment_options.watermelon.md5",
        reference_checksum = "config_checksums/config-references.watermelon.md5",
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
        ALIGNMENT_DIR + "/04-tophat/log/{sample}_tophat.log"
    threads: 8
    shell: 
            " module load watermelon_rnaseq && "
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
        ALIGNMENT_DIR + "/05-fastqc_align/log/{sample}_fastqc_tophat_align.log"
    shell:
        " module load watermelon_rnaseq && "
        "fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log} "

rule align_qc_metrics:
    input:
        sample_checksum = "config_checksums/config-samples.watermelon.md5",
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
        reference_checksum = "config_checksums/config-references.watermelon.md5",
        bams = ALIGNMENT_DIR + "/04-tophat/{sample}/{sample}_accepted_hits.bam",
        gtf = "references/gtf"
    output:
        DIFFEX_DIR + "/deseq2/01-htseq/{sample}_counts.txt"
    params:
        strand = rnaseq_snakefile_helper.check_strand_option("htseq", config["alignment_options"]["library_type"])
    log:
        DIFFEX_DIR + "/deseq2/01-htseq/log/{sample}_htseq_per_sample.log"
    shell:
        " module load watermelon_rnaseq &&"
        " python -m HTSeq.scripts.count "
        " -f bam "
        " -s {params.strand} "
        " -m "
        " intersection-nonempty "
        " -q {input.bams} "
        " {input.gtf} "
        " > {output} "
        " 2>&1 | tee {log}"

rule htseq_merge:
    input:
        sample_checksum = "config_checksums/config-samples.watermelon.md5",
        sample_count_files = expand("{diffex_dir}/deseq2/01-htseq/{sample}_counts.txt",
                                    diffex_dir=DIFFEX_DIR, sample=config["samples"])
    output:
        DIFFEX_DIR + "/deseq2/01-htseq/HTSeq_counts.txt"
    params:
        output_dir = DIFFEX_DIR + "/deseq2/01-htseq",
        input_dir = DIFFEX_DIR + "/deseq2/01-htseq"
    log:
        DIFFEX_DIR + "/deseq2/01-htseq/log/htseq_merge.log"
    shell:
       " perl {WATERMELON_SCRIPTS_DIR}/mergeHTSeqCountFiles.pl {params.input_dir} 2>&1 | tee {log}"

rule cuffdiff:
    input:
        sample_checksum = "config_checksums/phenotype_samples-{pheno}.watermelon.md5",
        comparison_checksum = "config_checksums/phenotype_comparisons-{pheno}.watermelon.md5",
        reference_checksum = "config_checksums/config-references.watermelon.md5",
        fasta_file = "references/bowtie2_index/genome.fa",
        gtf_file = "references/gtf",
        bam_files = expand("{alignment_dir}/04-tophat/{sample}/{sample}_accepted_hits.bam",
                            alignment_dir= ALIGNMENT_DIR,
                            sample=config["samples"])
    output:
        DIFFEX_DIR + "/08-cuffdiff/{pheno}/gene_exp.diff",
        DIFFEX_DIR + "/08-cuffdiff/{pheno}/isoform_exp.diff",
        DIFFEX_DIR + "/08-cuffdiff/{pheno}/read_groups.info"
    params:
        output_dir = DIFFEX_DIR + "/08-cuffdiff/{pheno}/",
        output_temp_dir = DIFFEX_DIR + "/08-cuffdiff/tmp_{pheno}/",
        labels = lambda wildcards : phenotypeManager.concatenated_comparison_values(',')[wildcards.pheno],
        samples = lambda wildcards : phenotypeManager.cuffdiff_samples(wildcards.pheno,
                                                                       ALIGNMENT_DIR + "/04-tophat/{sample_placeholder}/{sample_placeholder}_accepted_hits.bam"),
        strand = rnaseq_snakefile_helper.check_strand_option("tuxedo", config["alignment_options"]["library_type"])
    threads: 8
    log:
        DIFFEX_DIR + "/08-cuffdiff/log/{pheno}_cuffdiff.log"
    shell:
        " module load watermelon_rnaseq && "
        " cuffdiff -q "
        " -p {threads} "
        " -L {params.labels} "
        " --max-bundle-frags 999999999 "
        " --library-type {params.strand} "
        " -o {params.output_temp_dir} "
        " -b {input.fasta_file} "
        " -u -N "
        " --compatible-hits-norm "
        " {input.gtf_file} "
        " {params.samples} "
        " 2>&1 | tee {log} && "
        
        " mv {params.output_temp_dir}/* {params.output_dir} "

rule diffex_flip:
    input:
        gene_cuffdiff = DIFFEX_DIR + "/08-cuffdiff/{pheno}/gene_exp.diff",
        isoform_cuffdiff = DIFFEX_DIR + "/08-cuffdiff/{pheno}/isoform_exp.diff"
    output:
        gene_flip = DIFFEX_DIR + "/09-diffex_flip/{pheno}/gene_exp.flip.diff",
        isoform_flip = DIFFEX_DIR + "/09-diffex_flip/{pheno}/isoform_exp.flip.diff"
    params:
        comparisons = lambda wildcards: phenotypeManager.separated_comparisons(',')[wildcards.pheno]
    log: 
        DIFFEX_DIR + "/09-diffex_flip/log/{pheno}_diffex_flip.log"
    shell:
        " module purge && module load python/3.4.3 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_flip.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " {input.gene_cuffdiff} "
        " {output.gene_flip} "
        " {params.comparisons} "
        " 2>&1 | tee {log} && "
        
        "python {WATERMELON_SCRIPTS_DIR}/diffex_flip.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " {input.isoform_cuffdiff} "
        " {output.isoform_flip} "
        " {params.comparisons} "
        " 2>&1 | tee >>{log} "

rule diffex_flag:
    input:
        fold_change_checksum = "config_checksums/config-fold_change.watermelon.md5",
        cuffdiff_gene_exp = DIFFEX_DIR + "/09-diffex_flip/{pheno}/gene_exp.flip.diff",
        cuffdiff_isoform_exp = DIFFEX_DIR + "/09-diffex_flip/{pheno}/isoform_exp.flip.diff"
    output:
        gene_flagged = DIFFEX_DIR + "/10-diffex_flag/{pheno}/{pheno}_gene.flagged.txt",
        isoform_flagged = DIFFEX_DIR + "/10-diffex_flag/{pheno}/{pheno}_isoform.flagged.txt",
    params:
        fold_change = config["fold_change"]
    log:
        DIFFEX_DIR + "/10-diffex_flag/log/{pheno}_diffex_flag.log"
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
        genome_checksum = "config_checksums/config-genome.watermelon.md5",
        reference_checksum = "config_checksums/config-references.watermelon.md5",
        gene_diff_exp = DIFFEX_DIR + "/10-diffex_flag/{pheno}/{pheno}_gene.flagged.txt",
        isoform_diff_exp = DIFFEX_DIR + "/10-diffex_flag/{pheno}/{pheno}_isoform.flagged.txt",
        entrez_gene_info = "references/entrez_gene_info"
    output:
        gene_annot = DIFFEX_DIR + "/11-annotate_diffex_flag/{pheno}/{pheno}_gene.flagged.annot.txt",
        isoform_annot = DIFFEX_DIR + "/11-annotate_diffex_flag/{pheno}/{pheno}_isoform.flagged.annot.txt"
    params:
        output_dir = DIFFEX_DIR + "/11-annotate_diffex_flag/{pheno}",
        genome = config["genome"]
    log:
        DIFFEX_DIR + "/11-annotate_diffex_flag/log/{pheno}_annotate_diffex_flag.log"
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
        DIFFEX_DIR + "/08-cuffdiff/{pheno}/read_groups.info"
    output:
        DIFFEX_DIR + "/12-group_replicates/{pheno}/group_replicates.txt"
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
        genome_checksum = "config_checksums/config-genome.watermelon.md5",
        reference_checksum = "config_checksums/config-references.watermelon.md5",
        group_replicates = DIFFEX_DIR + "/12-group_replicates/{pheno}/group_replicates.txt",
        gtf_file = "references/gtf"
    output:
        DIFFEX_DIR + "/13-cummerbund/{pheno}/Plots",
        DIFFEX_DIR + "/13-cummerbund/{pheno}/Plots/{pheno}_boxplot.pdf",
        DIFFEX_DIR + "/13-cummerbund/{pheno}/{pheno}_repRawCounts.txt"
    params:
        cuff_diff_dir = DIFFEX_DIR + "/08-cuffdiff/{pheno}",
        output_dir = DIFFEX_DIR + "/13-cummerbund/{pheno}",
        genome = config["genome"],
    log:
         DIFFEX_DIR + "/13-cummerbund/log/{pheno}_cummerbund.log"
    shell:
        "module load watermelon_rnaseq && "
        "mkdir -p {params.output_dir}/Plots && "
        "Rscript {WATERMELON_SCRIPTS_DIR}/Run_cummeRbund.R "
        " baseDir={params.output_dir} "
        " cuffDiffDir={params.cuff_diff_dir} "
        " grpRepFile={input.group_replicates} "
        " gtfFile={input.gtf_file} "
        " genome={params.genome} "
        " 2>&1 | tee {log} && "
        " touch {params.output_dir}/Plots "
 
rule diffex_split:
    input:
        gene = "{diffex_dir}/11-annotate_diffex_flag/{phenotype_name}/{phenotype_name}_gene.flagged.annot.txt",
        isoform = "{diffex_dir}/11-annotate_diffex_flag/{phenotype_name}/{phenotype_name}_isoform.flagged.annot.txt",
    output:
        "{diffex_dir}/14-diffex_split/{phenotype_name}/{comparison}_gene.txt",
        "{diffex_dir}/14-diffex_split/{phenotype_name}/{comparison}_isoform.txt"
    params:
        output_dir = DIFFEX_DIR + "/14-diffex_split/{phenotype_name}",
        user_specified_comparison_list = lambda wildcards: phenotypeManager.separated_comparisons(',')[wildcards.phenotype_name],
    log:
        DIFFEX_DIR + "/14-diffex_split/log/{phenotype_name}_diffex_split.log"
    shell:
        "module purge && module load python/3.4.3 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_split.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " -o _gene.txt "
        " {input.gene} "
        " {params.output_dir} "
        " {params.user_specified_comparison_list} "
        " 2>&1 | tee {log} && "

        "module purge && module load python/3.4.3 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_split.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " -o _isoform.txt "
        " {input.isoform} "
        " {params.output_dir} "
        " {params.user_specified_comparison_list} "
        " 2>&1 | tee >>{log} "


rule last_split:
    input: 
        expand("{diffex_dir}/14-diffex_split/{phenotype_name}/{comparison}_gene.txt",
                zip,
                diffex_dir=repeat(DIFFEX_DIR, len(PHENOTYPE_NAMES)),
                phenotype_name=PHENOTYPE_NAMES,
                comparison=COMPARISON_GROUPS)
    output:
        touch(DIFFEX_DIR + "/14-diffex_split/last_split")


rule build_run_info:
    input: 
        DIFFEX_DIR + "/14-diffex_split/last_split",
        glossary = WATERMELON_SCRIPTS_DIR + "/glossary.txt"
        
    output: 
        DIFFEX_DIR + "/15-run_info/run_info.txt",
        DIFFEX_DIR + "/15-run_info/glossary.txt"

    run:
        command = 'module load watermelon_rnaseq; module list -t 2> {}'.format(output[0])
        subprocess.call(command, shell=True)
        with open(output[0], 'a') as run_info_file: 
            print('\n\nConfig\n', file=run_info_file)
            print(yaml.dump(config, default_flow_style=False), file=run_info_file)
        
        copyfile(WATERMELON_SCRIPTS_DIR + '/glossary.txt', DIFFEX_DIR + '/15-run_info/glossary.txt')
 
rule diffex_excel:
    input:
        gene = DIFFEX_DIR + "/14-diffex_split/{phenotype_name}/{comparison}_gene.txt",
        isoform = DIFFEX_DIR + "/14-diffex_split/{phenotype_name}/{comparison}_isoform.txt",
        glossary = DIFFEX_DIR + "/15-run_info/glossary.txt",
        run_info = DIFFEX_DIR + "/15-run_info/run_info.txt"
    output: 
        DIFFEX_DIR + "/16-diffex_excel/{phenotype_name}/{comparison}.xlsx"
    log:
        DIFFEX_DIR + "/16-diffex_excel/log/{phenotype_name}_diffex_excel.log"
    shell:
        "module purge && module load python/3.4.3 && "
        " python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py "
        " -g {input.gene}"
        " -i {input.isoform}"
        " --glossary {input.glossary} "
        " --info_filepath {input.run_info} "
        " {output} "
        " 2>&1 | tee {log} "

rule alignment_deliverables:
    input:
        raw_fastqc = expand("{alignment_dir}/03-fastqc_reads/{sample}_trimmed_R1_fastqc.html",
                                 alignment_dir=ALIGNMENT_DIR, sample=config["samples"]),
        align_fastqc = expand("{alignment_dir}/05-fastqc_align/{sample}_accepted_hits_fastqc.html",
                                 alignment_dir=ALIGNMENT_DIR, sample=config["samples"]),
        alignment_stats = ALIGNMENT_DIR + "/06-qc_metrics/alignment_stats.txt",
    output:
        raw_fastqc = expand(DIFFEX_DIR + "/Deliverables/alignment_deliverables/raw_reads_fastqc/{sample}_trimmed_R1_fastqc.html",
                            sample=config["samples"]),
        align_fastqc = expand(DIFFEX_DIR + "/Deliverables/alignment_deliverables/aligned_reads_fastqc/{sample}_accepted_hits_fastqc.html",
                            sample=config["samples"]),
        alignment_stats = DIFFEX_DIR + "/Deliverables/alignment_deliverables/alignment_stats.txt",
    params:
        raw_fastqc_input_dir    =  ALIGNMENT_DIR + "/03-fastqc_reads",
        raw_fastqc_output_dir   =  DIFFEX_DIR + "/Deliverables/alignment_deliverables/raw_reads_fastqc",
        align_fastqc_input_dir  =  ALIGNMENT_DIR + "/05-fastqc_align",
        align_fastqc_output_dir =  DIFFEX_DIR + "/Deliverables/alignment_deliverables/aligned_reads_fastqc",
    shell:
        " cp -r {params.raw_fastqc_input_dir}/* {params.raw_fastqc_output_dir} && "
        " cp -r {params.align_fastqc_input_dir}/* {params.align_fastqc_output_dir} && "
        " cp -r {input.alignment_stats} {output.alignment_stats} "

rule cuffdiff_deliverables:
    input:
        diffex_excel = expand("{diffex_dir}/16-diffex_excel/{phenotype_name}/{comparison}.xlsx", 
                                    zip,
                                    diffex_dir=repeat(DIFFEX_DIR, len(PHENOTYPE_NAMES)),
                                    phenotype_name=PHENOTYPE_NAMES, comparison=COMPARISON_GROUPS),
    output:
        diffex_excel = expand("{diffex_dir}/Deliverables/diffex_deliverables/cuffdiff_results/{phenotype_name}/{comparison}.xlsx", 
                                    zip,
                                    diffex_dir=repeat(DIFFEX_DIR, len(PHENOTYPE_NAMES)),
                                    phenotype_name=PHENOTYPE_NAMES, comparison=COMPARISON_GROUPS),
    params:
        diffex_excel_input_dir  =  DIFFEX_DIR +  "/16-diffex_excel",
        diffex_excel_output_dir =  DIFFEX_DIR + "/Deliverables/diffex_deliverables/cuffdiff_results"
    shell:
        " cp -r {params.diffex_excel_input_dir}/* {params.diffex_excel_output_dir} "

rule cummerbund_deliverables:
    input:
        diffex_raw_counts = "{diffex_dir}/13-cummerbund/{phenotype_name}/{phenotype_name}_repRawCounts.txt",
        plots = "{diffex_dir}/13-cummerbund/{phenotype_name}/Plots",
    output:
        diffex_raw_counts = "{diffex_dir}/Deliverables/diffex_deliverables/cummeRbund_results/{phenotype_name}_repRawCounts.txt",
        plots = "{diffex_dir}/Deliverables/diffex_deliverables/cummeRbund_results/{phenotype_name}_cummeRbund_plots"
    params:
        diffex_dir = "{diffex_dir}/Deliverables/diffex_deliverables",
    shell:
        " cp -r {input.diffex_raw_counts} {output.diffex_raw_counts} && "
        " find {output.diffex_raw_counts} | xargs -I ^ touch ^ && "
        " cp -r {input.plots}  {output.plots} && "
        " find {output.plots} | xargs -I ^ touch ^ "

rule deseq2_metadata_contrasts:
    input:
        sample_checksum      = "config_checksums/config-samples.watermelon.md5",
        comparison_checksum  = "config_checksums/config-comparisons.watermelon.md5",
        phenotype_checksum   = "config_checksums/config-phenotypes.watermelon.md5",
        main_factor_checksum = "config_checksums/config-main_factors.watermelon.md5"
    output:
        sample_metadata = DIFFEX_DIR + "/deseq2/02-metadata_contrasts/sample_metadata.txt",
        contrasts = DIFFEX_DIR + "/deseq2/02-metadata_contrasts/contrasts.txt"
    params:
        output_dir =DIFFEX_DIR + "/deseq2/03-deseq2_diffex"
    run:
        deseq2_helper.build_sample_metadata(config, output.sample_metadata)
        deseq2_helper.build_contrasts(config, params.output_dir, output.contrasts)

        
