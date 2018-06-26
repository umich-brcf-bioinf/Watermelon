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

WATERMELON_CONFIG_DIR = os.path.join(os.environ.get('WATERMELON_CONFIG_DIR', 'config'), '')
WATERMELON_SCRIPTS_DIR = os.path.join(os.environ.get('WATERMELON_SCRIPTS_DIR', 'scripts'), '')

INPUT_DIR = os.path.join(config.get("input_dir", "inputs"), "")
ALIGNMENT_DIR = os.path.join(config.get("alignment_output_dir", "alignment_results"), "")
DIFFEX_DIR = os.path.join(config.get("diffex_output_dir", "diffex_results"), "")
DELIVERABLES_DIR = os.path.join(config.get("deliverables_output_dir", "deliverables"), "")
DESEQ2_DIR = os.path.join(DIFFEX_DIR, "deseq2", "")
TUXEDO_DIR = os.path.join(DIFFEX_DIR, "tuxedo", "")
CONFIG_CHECKSUMS_DIR = os.path.join(".config_checksums", "")

COMPARISON_INFIX = '_v_'
DELIMITER = '^'
SAMPLES_KEY = 'samples'
PHENOTYPES_KEY = 'phenotypes'
COMPARISONS_KEY ='comparisons'

phenotypeManager = rnaseq_snakefile_helper.PhenotypeManager(config,
                                                            DELIMITER,
                                                            COMPARISON_INFIX)

ALL_PHENOTYPE_NAMES = phenotypeManager.phenotypes_comparisons_all_tuple.phenotypes
ALL_COMPARISON_GROUPS = phenotypeManager.phenotypes_comparisons_all_tuple.comparisons

REPLICATE_PHENOTYPE_NAMES = phenotypeManager.phenotypes_comparisons_replicates_tuple.phenotypes
REPLICATE_COMPARISON_GROUPS = phenotypeManager.phenotypes_comparisons_replicates_tuple.comparisons


rnaseq_snakefile_helper.init_references(config["references"])
rnaseq_snakefile_helper.checksum_reset_all(CONFIG_CHECKSUMS_DIR,
                                           config=config,
                                           phenotype_comparisons=config['comparisons'],
                                           phenotype_samples=phenotypeManager.phenotype_sample_list)

SAMPLE_READS = rnaseq_snakefile_helper.flattened_sample_reads(config['input_dir'], config[SAMPLES_KEY])

if 'fastq_screen' in config:
    FASTQ_SCREEN_CONFIG = config['fastq_screen']
    FASTQ_SCREEN_ALIGNMENT = [
            rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "03-fastq_screen/multi_species/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
            rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "03-fastq_screen/biotype/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS)
        ]
    FASTQ_SCREEN_DELIVERABLES = [
            rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/fastq_screen/multi_species/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
            rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/fastq_screen/biotype/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
        ]
else:
    FASTQ_SCREEN_CONFIG = defaultdict(str)
    FASTQ_SCREEN_ALIGNMENT = []
    FASTQ_SCREEN_DELIVERABLES = []

DESEQ2_ALL = []
if REPLICATE_PHENOTYPE_NAMES:
    DESEQ2_ALL = [
        expand(DESEQ2_DIR + "01-htseq/{sample}_counts.txt",
               sample=config[SAMPLES_KEY]),
        DESEQ2_DIR + "01-htseq/htseq_merged.txt",
        DESEQ2_DIR + "02-metadata_contrasts/sample_metadata.txt",
        DESEQ2_DIR + "02-metadata_contrasts/contrasts.txt",
        expand(DESEQ2_DIR + "03-deseq2_diffex/gene_lists/{phenotype_name}/{comparison}.txt",
               zip,
               phenotype_name=REPLICATE_PHENOTYPE_NAMES,
               comparison=REPLICATE_COMPARISON_GROUPS),
        expand(DESEQ2_DIR + "04-annotation/{phenotype_name}/{comparison}.annot.txt",
               zip,
               phenotype_name=REPLICATE_PHENOTYPE_NAMES,
               comparison=REPLICATE_COMPARISON_GROUPS),
        expand(DESEQ2_DIR + "06-excel/{phenotype_name}/{comparison}.xlsx",
                zip,
                phenotype_name=REPLICATE_PHENOTYPE_NAMES,
                comparison=REPLICATE_COMPARISON_GROUPS),
        DESEQ2_DIR + "07-summary/deseq2_summary.txt",
        DESEQ2_DIR + "07-summary/deseq2_summary.xlsx",
        DELIVERABLES_DIR + "deseq2",
        ]

OPTIONAL_ALL = DESEQ2_ALL + FASTQ_SCREEN_ALIGNMENT + FASTQ_SCREEN_DELIVERABLES

rule all:
    input:
        rnaseq_snakefile_helper.expand_sample_read_endedness(\
            ALIGNMENT_DIR + "03-fastqc_reads/{sample}_trimmed_{read_endedness}_fastqc.html",
            SAMPLE_READS),
        expand(ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_accepted_hits.bam",
                sample=config[SAMPLES_KEY]),
        expand(ALIGNMENT_DIR + "05-fastqc_align/{sample}_accepted_hits_fastqc.html",
                sample=config[SAMPLES_KEY]),
        ALIGNMENT_DIR + "06-qc/alignment_qc.html",
        rnaseq_snakefile_helper.expand_sample_read_endedness(\
            DELIVERABLES_DIR + "alignment/sequence_reads_fastqc/{sample}_trimmed_{read_endedness}_fastqc.html",
            SAMPLE_READS),
        expand(DELIVERABLES_DIR + "alignment/aligned_reads_fastqc/{sample}_accepted_hits_fastqc.html",
                sample=config[SAMPLES_KEY]),
        DELIVERABLES_DIR + "alignment/alignment_qc.html",
        expand(TUXEDO_DIR + "01-cuffdiff/{phenotype}/gene_exp.diff",
               phenotype=sorted(config[COMPARISONS_KEY].keys())),
        expand(TUXEDO_DIR + "03-flag/{phenotype}/{phenotype}_gene.flagged.txt",
               phenotype=sorted(config[COMPARISONS_KEY].keys())),
        expand(TUXEDO_DIR + "03-flag/{phenotype}/{phenotype}_isoform.flagged.txt",
                phenotype=sorted(config[COMPARISONS_KEY].keys())),
        expand(TUXEDO_DIR + "04-annotate/{phenotype}/{phenotype}_gene.flagged.annot.txt",
               phenotype=sorted(config[COMPARISONS_KEY].keys())),
        expand(TUXEDO_DIR + "04-annotate/{phenotype}/{phenotype}_isoform.flagged.annot.txt",
               phenotype=sorted(config[COMPARISONS_KEY].keys())),
        expand(TUXEDO_DIR + "06-cummerbund/{pheno}/Plots/{pheno}_boxplot.pdf",
               pheno=sorted(config[COMPARISONS_KEY].keys())),
        expand(TUXEDO_DIR + "06-cummerbund/{pheno}/{pheno}_repRawCounts.txt",
               pheno=sorted(config[COMPARISONS_KEY].keys())),
        expand(TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_gene.txt",
                zip,
                phenotype_name=ALL_PHENOTYPE_NAMES,
                comparison=ALL_COMPARISON_GROUPS),
        expand(TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_isoform.txt",
                zip,
                phenotype_name=ALL_PHENOTYPE_NAMES,
                comparison=ALL_COMPARISON_GROUPS),
        TUXEDO_DIR + "07-split/last_split",
        expand(TUXEDO_DIR + "09-excel/{phenotype_name}/{comparison}.xlsx",
                zip,
                phenotype_name=ALL_PHENOTYPE_NAMES,
                comparison=ALL_COMPARISON_GROUPS),
        TUXEDO_DIR + "10-summary/tuxedo_summary.txt",
        TUXEDO_DIR + "10-summary/tuxedo_summary.xlsx",
        DELIVERABLES_DIR + "tuxedo",
        DELIVERABLES_DIR + "combined_summary.xlsx",
        *OPTIONAL_ALL

include: 'rules/align_concat_reads.smk'
include: 'rules/align_cutadapt_SE.smk'
include: 'rules/align_cutadapt_PE.smk'
include: 'rules/align_fastq_screen_biotype.smk'
include: 'rules/align_fastq_screen_multi_species.smk'

rule align_fastqc_trimmed_reads:
    input:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_{read_endedness}.fastq.gz"
    output:
        touch(ALIGNMENT_DIR + "03-fastqc_reads/reads_fastq.done"),
        ALIGNMENT_DIR + "03-fastqc_reads/{sample}_trimmed_{read_endedness}_fastqc.html"
    log:
        ALIGNMENT_DIR + "03-fastqc_reads/.log/{sample}_trimmed_{read_endedness}_fastqc.log"
    params:
        fastqc_dir = ALIGNMENT_DIR + "03-fastqc_reads"
    shell:
        "module purge && module load watermelon_dependencies && "
        "fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log}"

rule align_insert_size_PE:
    input:
        fasta_file = "references/bowtie2_index/genome.fa",
        bbmap_ref = "references/bbmap_ref",
        raw_fastq_R1 = ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_R1_PE.fastq.gz",
        raw_fastq_R2 = ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_R2_PE.fastq.gz",
    output:
        ALIGNMENT_DIR + "04-insert_size/{sample}_read_stats.txt"
    params:
        sample_id = "{sample}",
        downsample_readcount = config["insert_size"]["downsample_readcount"],
        pairlen = config["insert_size"]["pair_length"],
        temp_dir = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}",
        downsampled_fastq_r1 = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}/{sample}_R1.downsampled.fastq.gz",
        downsampled_fastq_r2 = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}/{sample}_R2.downsampled.fastq.gz",
        ihist_file = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}/{sample}_ihist.txt",
        lhist_file = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}/{sample}_lhist.txt",
        mapped_file = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}/{sample}.sam",
        output_dir = ALIGNMENT_DIR + "04-insert_size/",
    threads:
        8
    resources:
        memoryInGb=32
    log:
        ALIGNMENT_DIR + "04-insert_size/.log/{sample}_read_stats.log"
    shell:
        '''(module purge && module load watermelon_dependencies &&
        mkdir -p {params.temp_dir} &&
        JAVA_OPTS="-XX:ParallelGCThreads=2" &&
        reformat.sh t={threads} \
            sampleseed=42 \
            samplereadstarget={params.downsample_readcount} \
            in1={input.raw_fastq_R1} in2={input.raw_fastq_R2} \
            out={params.downsampled_fastq_r1} out2={params.downsampled_fastq_r2} &&
        bbmap.sh -Xmx{resources.memoryInGb}g t={threads} \
            semiperfectmode=t \
            requirecorrectstrand=f \
            pairlen={params.pairlen} \
            path={input.bbmap_ref} \
            in1={params.downsampled_fastq_r1} \
            in2={params.downsampled_fastq_r2} \
            ihist={params.ihist_file} \
            lhist={params.lhist_file} \
            out={params.mapped_file} &&
        module purge && module load python/3.6.1 &&
        python {WATERMELON_SCRIPTS_DIR}/read_stats_bbmap_single.py \
            --input_dir {params.temp_dir} \
            --sample_id {params.sample_id} \
            --output_dir {params.output_dir}
        ) 2>&1 | tee {log}'''

rule align_create_transcriptome_index:
    input:
        alignment_options_checksum =  CONFIG_CHECKSUMS_DIR + "config-alignment_options.watermelon.md5",
        reference_checksum =  CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        gtf = "references/gtf",
        bowtie2_index_dir = "references/bowtie2_index"
    output:
        ALIGNMENT_DIR + "04-tophat/transcriptome_index/transcriptome.fa"
    params:
        transcriptome_dir = "transcriptome_index",
        temp_dir =  ALIGNMENT_DIR + "04-tophat/.tmp",
        output_dir = ALIGNMENT_DIR + "04-tophat",
        strand = rnaseq_snakefile_helper.strand_option_tophat(config["alignment_options"]["library_type"])
    log:
        ALIGNMENT_DIR + "04-tophat/.log/create_transcriptome_index.log"
    shell:
        '''(module purge && module load watermelon_dependencies &&
        mkdir -p {params.temp_dir} &&
        rm -rf {params.temp_dir}/* &&
        tophat -G {input.gtf} \
            --library-type {params.strand} \
            --transcriptome-index={params.temp_dir}/transcriptome_index/transcriptome \
            {input.bowtie2_index_dir}/genome &&
        rm -rf {params.output_dir}/{params.transcriptome_dir} &&
        mv {params.temp_dir}/{params.transcriptome_dir} {params.output_dir} &&
        mv tophat_out {params.output_dir}/{params.transcriptome_dir}/ &&
        touch {params.output_dir}/{params.transcriptome_dir}/*
        ) 2>&1 | tee {log}'''

rule align_build_tophat_sample_options:
    input:
        read_stats = lambda wildcards: rnaseq_snakefile_helper.expand_read_stats_if_paired(\
                ALIGNMENT_DIR + "04-insert_size/{sample}_read_stats.txt",
                SAMPLE_READS,
                wildcards.sample)
    output:
        tophat_sample_options = ALIGNMENT_DIR + "04-tophat/{sample}/{sample}.tophat_sample_options",
    run:
        options = ''
        with open(str(output['tophat_sample_options']), 'w') as output:
            if 'read_stats' in input:
                options = rnaseq_snakefile_helper.tophat_paired_end_flags(str(input['read_stats']))
            output.write(options)

rule align_tophat:
    input:
        alignment_options_checksum = CONFIG_CHECKSUMS_DIR + "config-alignment_options.watermelon.md5",
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        transcriptome_fasta = ALIGNMENT_DIR + "04-tophat/transcriptome_index/transcriptome.fa",
        bowtie2_index_dir = "references/bowtie2_index",
        fastq_files = lambda wildcards: rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_{read_endedness}.fastq.gz",
                SAMPLE_READS,
                wildcards.sample),
        tophat_sample_options = ALIGNMENT_DIR + "04-tophat/{sample}/{sample}.tophat_sample_options",
    output:
        ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_accepted_hits.bam",
        ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_align_summary.txt",
    params:
        transcriptome_index = ALIGNMENT_DIR + "04-tophat/transcriptome_index/transcriptome",
        tophat_alignment_options = rnaseq_snakefile_helper.tophat_options(config["alignment_options"]),
        strand = rnaseq_snakefile_helper.strand_option_tophat(config["alignment_options"]["library_type"]),
        tophat_dir = ALIGNMENT_DIR + "04-tophat",
        sample = lambda wildcards: wildcards.sample,
        paired_end_flags = lambda wildcards: rnaseq_snakefile_helper.tophat_paired_end_flags(\
                ALIGNMENT_DIR + "04-insert_size/{sample}_read_stats.txt".format(sample=wildcards.sample)),
    log:
        ALIGNMENT_DIR + "04-tophat/.log/{sample}_tophat.log"
    threads: 8
    shell:
        '''(module purge && module load watermelon_dependencies &&
        TOPHAT_SAMPLE_OPTIONS=$(<{input.tophat_sample_options}) &&
        set -x &&
        tophat -p {threads} \
             --b2-very-sensitive \
             --no-coverage-search \
             --library-type {params.strand} \
             -I 500000 \
             --transcriptome-index={params.transcriptome_index} \
             $TOPHAT_SAMPLE_OPTIONS \
             {params.tophat_alignment_options} \
             -o {params.tophat_dir}/{params.sample} \
             {input.bowtie2_index_dir}/genome \
             {input.fastq_files} && \
        mv {params.tophat_dir}/{params.sample}/accepted_hits.bam {params.tophat_dir}/{params.sample}/{params.sample}_accepted_hits.bam &&
        mv {params.tophat_dir}/{params.sample}/align_summary.txt {params.tophat_dir}/{params.sample}/{params.sample}_align_summary.txt
        ) 2>&1 | tee {log}'''

rule align_fastqc_tophat_align:
    input:
        ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_accepted_hits.bam"
    output:
        touch(ALIGNMENT_DIR + "05-fastqc_align/align_fastq.done"),
        ALIGNMENT_DIR + "05-fastqc_align/{sample}_accepted_hits_fastqc.html"
    params:
        fastqc_dir =  ALIGNMENT_DIR + "05-fastqc_align"
    log:
        ALIGNMENT_DIR + "05-fastqc_align/.log/{sample}_fastqc_tophat_align.log"
    shell:
        "module purge && module load watermelon_dependencies && "
        "fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log} "

rule align_qc:
    input:
        sample_checksum = CONFIG_CHECKSUMS_DIR + "config-samples.watermelon.md5",
        raw_read_fastq_files = rnaseq_snakefile_helper.expand_sample_read_endedness(\
            ALIGNMENT_DIR + "03-fastqc_reads/{sample}_trimmed_{read_endedness}_fastqc.html",
            SAMPLE_READS),
        align_summary_files = expand(ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_align_summary.txt",
                                     sample=config["samples"]),
        align_fastq_files = expand(ALIGNMENT_DIR + "05-fastqc_align/{sample}_accepted_hits_fastqc.html",
                                     sample=config["samples"]),
        fastq_screen_alignment = FASTQ_SCREEN_ALIGNMENT,
    output:
        ALIGNMENT_DIR + "06-qc/alignment_qc.html"
    params:
        output_dir = ALIGNMENT_DIR + "06-qc/",
        output_filename = "alignment_qc.html",
        multiqc_config_filename = WATERMELON_CONFIG_DIR + "multiqc_config.yaml",
    log:
        ALIGNMENT_DIR + "06-qc/.log/align_qc.log"
    shell:
        '''(module purge && module load watermelon_dependencies &&
        echo 'watermelon|version|multiqc|'`multiqc --version | cut -d' ' -f2-` &&
        multiqc --force \
            --exclude cutadapt \
            --exclude bowtie2 \
            --config {params.multiqc_config_filename} \
            --outdir {params.output_dir} \
            --filename {params.output_filename} \
            alignment_results
        ) 2>&1 | tee {log} '''

rule align_deliverables_alignment:
    input:
        rnaseq_snakefile_helper.expand_sample_read_endedness(
                ALIGNMENT_DIR + "03-fastqc_reads/{sample}_trimmed_{read_endedness}_fastqc.html",
                SAMPLE_READS),
        expand(ALIGNMENT_DIR + "05-fastqc_align/{sample}_accepted_hits_fastqc.html",
                sample=config["samples"]),
        alignment_stats = ALIGNMENT_DIR + "06-qc/alignment_qc.html",
    output:
        rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/sequence_reads_fastqc/{sample}_trimmed_{read_endedness}_fastqc.html",
                SAMPLE_READS),
        expand(DELIVERABLES_DIR + "alignment/aligned_reads_fastqc/{sample}_accepted_hits_fastqc.html",
                sample=config["samples"]),
        alignment_stats = DELIVERABLES_DIR + "alignment/alignment_qc.html",
    params:
        raw_fastqc_input_dir    =  ALIGNMENT_DIR + "03-fastqc_reads",
        raw_fastqc_output_dir   =  DELIVERABLES_DIR + "alignment/sequence_reads_fastqc",
        align_fastqc_input_dir  =  ALIGNMENT_DIR + "05-fastqc_align",
        align_fastqc_output_dir =  DELIVERABLES_DIR + "alignment/aligned_reads_fastqc",
    shell:
        "cp -r {params.raw_fastqc_input_dir}/* {params.raw_fastqc_output_dir} && "
        "cp -r {params.align_fastqc_input_dir}/* {params.align_fastqc_output_dir} && "
        "cp -r {input.alignment_stats} {output.alignment_stats} "

rule align_deliverables_fastq_screen:
    input:
        rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "03-fastq_screen/multi_species/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
        rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "03-fastq_screen/biotype/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
    output:
        rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/fastq_screen/multi_species/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
        rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/fastq_screen/biotype/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
    params:
        source_multi_species    =  ALIGNMENT_DIR + "03-fastq_screen/multi_species/",
        dest_multi_species    =  DELIVERABLES_DIR + "alignment/fastq_screen/multi_species/",
        source_biotype    =  ALIGNMENT_DIR + "03-fastq_screen/biotype/",
        dest_biotype    =  DELIVERABLES_DIR + "alignment/fastq_screen/biotype/",
    shell:
        "cp -r {params.source_multi_species}/*_screen.html {params.dest_multi_species} && "
        "cp -r {params.source_biotype}/*_screen.html {params.dest_biotype}"


rule tuxedo_cuffdiff:
    input:
        sample_checksum = CONFIG_CHECKSUMS_DIR + "phenotype_samples-{pheno}.watermelon.md5",
        comparison_checksum = CONFIG_CHECKSUMS_DIR + "phenotype_comparisons-{pheno}.watermelon.md5",
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        fasta_file = "references/bowtie2_index/genome.fa",
        gtf_file = "references/gtf",
        bam_files = expand(ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_accepted_hits.bam",
                           sample=config[SAMPLES_KEY])
    output:
        TUXEDO_DIR + "01-cuffdiff/{pheno}/gene_exp.diff",
        TUXEDO_DIR + "01-cuffdiff/{pheno}/isoform_exp.diff",
        TUXEDO_DIR + "01-cuffdiff/{pheno}/read_groups.info"
    params:
        output_dir = TUXEDO_DIR + "01-cuffdiff/{pheno}",
        labels = lambda wildcards : phenotypeManager.concatenated_comparison_values(',')[wildcards.pheno],
        samples = lambda wildcards : phenotypeManager.cuffdiff_samples(wildcards.pheno,
                                                                       ALIGNMENT_DIR + "04-tophat/{sample_placeholder}/{sample_placeholder}_accepted_hits.bam"),
        strand = rnaseq_snakefile_helper.strand_option_tophat(config["alignment_options"]["library_type"])
    threads: 8
    log:
        TUXEDO_DIR + "01-cuffdiff/.log/{pheno}_cuffdiff.log"
    shell:
        "rm -rf {params.output_dir} {params.output_dir}.tmp &&"
        "module purge && module load watermelon_dependencies && "
        "cuffdiff -q "
        " -p {threads} "
        " -L {params.labels} "
        " --max-bundle-frags 999999999 "
        " --library-type {params.strand} "
        " -o {params.output_dir}.tmp "
        " -b {input.fasta_file} "
        " -u -N "
        " --compatible-hits-norm "
        " {input.gtf_file} "
        " {params.samples} "
        " 2>&1 | tee {log} && "
        "mv {params.output_dir}.tmp {params.output_dir} "

rule tuxedo_flip:
    input:
        gene_cuffdiff = TUXEDO_DIR + "01-cuffdiff/{pheno}/gene_exp.diff",
        isoform_cuffdiff = TUXEDO_DIR + "01-cuffdiff/{pheno}/isoform_exp.diff"
    output:
        gene_flip = TUXEDO_DIR + "02-flip/{pheno}/gene_exp.flip.diff",
        isoform_flip = TUXEDO_DIR + "02-flip/{pheno}/isoform_exp.flip.diff"
    params:
        comparisons = lambda wildcards: phenotypeManager.separated_comparisons(',')[wildcards.pheno]
    log:
        TUXEDO_DIR + "02-flip/.log/{pheno}_flip.log"
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_flip.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " {input.gene_cuffdiff} "
        " {output.gene_flip} "
        " {params.comparisons} "
        " 2>&1 | tee {log} && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_flip.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " {input.isoform_cuffdiff} "
        " {output.isoform_flip} "
        " {params.comparisons} "
        " 2>&1 | tee >>{log} "

rule tuxedo_flag:
    input:
        fold_change_checksum = CONFIG_CHECKSUMS_DIR + "config-fold_change.watermelon.md5",
        cuffdiff_gene_exp = TUXEDO_DIR + "02-flip/{pheno}/gene_exp.flip.diff",
        cuffdiff_isoform_exp = TUXEDO_DIR + "02-flip/{pheno}/isoform_exp.flip.diff"
    output:
        gene_flagged = TUXEDO_DIR + "03-flag/{pheno}/{pheno}_gene.flagged.txt",
        isoform_flagged = TUXEDO_DIR + "03-flag/{pheno}/{pheno}_isoform.flagged.txt",
    params:
        fold_change = config["fold_change"]
    log:
        TUXEDO_DIR + "03-flag/.log/{pheno}_tuxedo_flag.log"
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_flag.py "
        " -f {params.fold_change} "
        " {input.cuffdiff_gene_exp} "
        " {output.gene_flagged} "
        " 2>&1 | tee {log} && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_flag.py "
        " -f {params.fold_change} "
        " {input.cuffdiff_isoform_exp} "
        " {output.isoform_flagged} "
        " 2>&1 | tee >>{log} "

rule tuxedo_annotate:
    input:
        genome_checksum = CONFIG_CHECKSUMS_DIR + "config-genome.watermelon.md5",
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        gene_diff_exp = TUXEDO_DIR + "03-flag/{pheno}/{pheno}_gene.flagged.txt",
        isoform_diff_exp = TUXEDO_DIR + "03-flag/{pheno}/{pheno}_isoform.flagged.txt",
        entrez_gene_info = "references/entrez_gene_info"
    output:
        gene_annot = TUXEDO_DIR + "04-annotate/{pheno}/{pheno}_gene.flagged.annot.txt",
        isoform_annot = TUXEDO_DIR + "04-annotate/{pheno}/{pheno}_isoform.flagged.annot.txt"
    params:
        output_dir = TUXEDO_DIR + "04-annotate/{pheno}",
        genome = config["genome"]
    log:
        TUXEDO_DIR + "04-annotate/.log/{pheno}_annotate.log"
    shell:
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_annotate.py "
        " -i {input.entrez_gene_info} "
        " -e {input.gene_diff_exp} "
        " -g {params.genome} "
        " -o {params.output_dir} "
        " 2>&1 | tee {log} && "

        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_annotate.py "
        " -i {input.entrez_gene_info} "
        " -e {input.isoform_diff_exp} "
        " -g {params.genome} "
        " -o {params.output_dir} "
        " 2>&1 | tee >>{log} "

rule tuxedo_group_replicates:
    input:
        TUXEDO_DIR + "01-cuffdiff/{pheno}/read_groups.info"
    output:
        TUXEDO_DIR + "05-group_replicates/{pheno}/group_replicates.txt"
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

rule tuxedo_cummerbund:
    input:
        genome_checksum = CONFIG_CHECKSUMS_DIR + "config-genome.watermelon.md5",
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        group_replicates = TUXEDO_DIR + "05-group_replicates/{pheno}/group_replicates.txt",
        gtf_file = "references/gtf"
    output:
        TUXEDO_DIR + "06-cummerbund/{pheno}/Plots",
        TUXEDO_DIR + "06-cummerbund/{pheno}/Plots/{pheno}_boxplot.pdf",
        TUXEDO_DIR + "06-cummerbund/{pheno}/{pheno}_repRawCounts.txt"
    params:
        cuff_diff_dir = TUXEDO_DIR + "01-cuffdiff/{pheno}",
        output_dir = TUXEDO_DIR + "06-cummerbund/{pheno}",
        genome = config["genome"],
    log:
         TUXEDO_DIR + "06-cummerbund/.log/{pheno}_cummerbund.log"
    shell:
        "module purge && module load watermelon_dependencies && "
        "mkdir -p {params.output_dir}/Plots && "
        "Rscript {WATERMELON_SCRIPTS_DIR}/Run_cummeRbund.R "
        " baseDir={params.output_dir} "
        " cuffDiffDir={params.cuff_diff_dir} "
        " grpRepFile={input.group_replicates} "
        " gtfFile={input.gtf_file} "
        " genome={params.genome} "
        " 2>&1 | tee {log} && "
        "touch {params.output_dir}/Plots "

rule tuxedo_split:
    input:
        gene = TUXEDO_DIR + "04-annotate/{phenotype_name}/{phenotype_name}_gene.flagged.annot.txt",
        isoform = TUXEDO_DIR + "04-annotate/{phenotype_name}/{phenotype_name}_isoform.flagged.annot.txt",
    output:
        TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_gene.txt",
        TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_isoform.txt",
    params:
        output_dir = TUXEDO_DIR + "07-split/{phenotype_name}",
        user_specified_comparison_list = lambda wildcards: phenotypeManager.separated_comparisons(',')[wildcards.phenotype_name],
    log:
        TUXEDO_DIR + "07-split/.log/{phenotype_name}_tuxedo_split.log"
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_split.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " -o _gene.txt "
        " {input.gene} "
        " {params.output_dir} "
        " {params.user_specified_comparison_list} "
        " 2>&1 | tee {log} && "

        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_split.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " -o _isoform.txt "
        " {input.isoform} "
        " {params.output_dir} "
        " {params.user_specified_comparison_list} "
        " 2>&1 | tee >>{log} "


rule tuxedo_last_split:
    input:
        expand(TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_gene.txt", \
               zip, \
               phenotype_name=ALL_PHENOTYPE_NAMES, \
               comparison=ALL_COMPARISON_GROUPS)
    output:
        touch(TUXEDO_DIR + "07-split/last_split")


rule tuxedo_run_info:
    input:
        TUXEDO_DIR + "07-split/last_split",
        glossary = WATERMELON_SCRIPTS_DIR + "tuxedo_glossary.txt"
    output:
        run_info=TUXEDO_DIR + "08-run_info/run_info.txt",
        glossary=TUXEDO_DIR + "08-run_info/glossary.txt"
    run:
        command = ('module load watermelon_dependencies && '
                   'module list -t 2> {}').format(output['run_info'])
        subprocess.call(command, shell=True)
        with open(output['run_info'], 'a') as run_info_file:
            print('\n\nConfig\n', file=run_info_file)
            print(yaml.dump(config, default_flow_style=False),
                  file=run_info_file)
        copyfile(input['glossary'], output['glossary'])

rule tuxedo_excel:
    input:
        gene = TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_gene.txt",
        isoform = TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_isoform.txt",
        glossary = TUXEDO_DIR + "08-run_info/glossary.txt",
        run_info = TUXEDO_DIR + "08-run_info/run_info.txt"
    output:
        TUXEDO_DIR + "09-excel/{phenotype_name}/{comparison}.xlsx"
    log:
        TUXEDO_DIR + "09-excel/.log/{phenotype_name}_diffex_excel.log"
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py "
        " -g {input.gene}"
        " -i {input.isoform}"
        " --glossary {input.glossary} "
        " --info_filepath {input.run_info} "
        " {output} "
        " 2>&1 | tee {log} "

rule tuxedo_summary:
    input:
        input_files = expand(TUXEDO_DIR + "07-split/{phenotype}/{comparison}_gene.txt",
                             zip,
                             phenotype=ALL_PHENOTYPE_NAMES,
                             comparison=ALL_COMPARISON_GROUPS) +
                       expand(TUXEDO_DIR + "07-split/{phenotype}/{comparison}_isoform.txt",
                              zip,
                              phenotype=ALL_PHENOTYPE_NAMES,
                              comparison=ALL_COMPARISON_GROUPS)
    output:
        summary_txt = TUXEDO_DIR + "10-summary/tuxedo_summary.txt",
        summary_xlsx = TUXEDO_DIR + "10-summary/tuxedo_summary.xlsx",
    log:
        TUXEDO_DIR + "10-summary/.log/tuxedo_summary.log"
    params:
        output_dir = TUXEDO_DIR + "10-summary/",
    shell:
        '''(module purge && module load python/3.6.1 &&
        set -v &&
        python {WATERMELON_SCRIPTS_DIR}/diffex_summary.py \
            --annotation_column gene_id \
            --annotation_null . \
            --diffex_call_column diff_exp \
            --diffex_call_pass Yes \
            --trim_suffix .txt \
            --output_file {output.summary_txt} \
            --output_xlsx {output.summary_xlsx} \
            {input.input_files} &&
        touch {params.output_dir}
        ) 2>&1 | tee {log} '''

rule deliverables_tuxedo:
    input:
        excel = expand(TUXEDO_DIR + "09-excel/{phenotype_name}/{comparison}.xlsx",
                       zip,
                       phenotype_name=ALL_PHENOTYPE_NAMES,
                       comparison=ALL_COMPARISON_GROUPS),
        diffex_raw_counts = expand(TUXEDO_DIR + "06-cummerbund/{phenotype_name}/{phenotype_name}_repRawCounts.txt",
                                   phenotype_name=ALL_PHENOTYPE_NAMES),
        plots = expand(TUXEDO_DIR + "06-cummerbund/{phenotype_name}/Plots",
                       phenotype_name=ALL_PHENOTYPE_NAMES),
        summary_txt = TUXEDO_DIR + "10-summary/tuxedo_summary.txt",
        summary_xlsx = TUXEDO_DIR + "10-summary/tuxedo_summary.xlsx",
    output:
        deliverables_dir = DELIVERABLES_DIR + "tuxedo",
        summary_gene_lists = DELIVERABLES_DIR + "tuxedo/gene_lists/tuxedo_summary.txt",
    params:
        tmp_dir = DELIVERABLES_DIR + "tuxedo.tmp",
        source_gene_list_dir = TUXEDO_DIR +  "/09-excel",
        source_counts_dir =  TUXEDO_DIR +  "/06-cummerbund",
        source_plots_dir =  TUXEDO_DIR +  "/06-cummerbund",
    shell:
        """rm -rf {output.deliverables_dir} {params.tmp_dir}
        mkdir -p {params.tmp_dir}/counts {params.tmp_dir}/plots
        cp `find {params.source_counts_dir} -maxdepth 2 -name '*_repRawCounts.txt'` {params.tmp_dir}/counts
        (for phenotype in `ls {params.source_plots_dir}`; do cp -r {params.source_plots_dir}/${{phenotype}}/Plots {params.tmp_dir}/plots/$phenotype; done)
        cp -r {params.source_gene_list_dir} {params.tmp_dir}/gene_lists
        cp {input.summary_txt} {params.tmp_dir}/gene_lists
        cp {input.summary_xlsx} {params.tmp_dir}/gene_lists
        mv {params.tmp_dir} {output.deliverables_dir} """

def gene_summaries():
    summaries = [DELIVERABLES_DIR + "tuxedo/gene_lists/tuxedo_summary.txt"]
    if DESEQ2_ALL:
        summaries.append(DELIVERABLES_DIR + "deseq2/gene_lists/deseq2_summary.txt")
    return summaries

rule deliverables_combined_summary:
    input:
        lambda x: gene_summaries()
    output:
        txt = DELIVERABLES_DIR + "combined_summary.txt",
        xlsx = DELIVERABLES_DIR + "combined_summary.xlsx",
    params:
        base_name = DELIVERABLES_DIR + "combined_summary"
    shell:
        '''python {WATERMELON_SCRIPTS_DIR}/combine_summaries.py \
            --output_base {params.base_name} {input}
        '''

rule deseq2_htseq:
    input:
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        bams = ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_accepted_hits.bam",
        gtf = "references/gtf"
    output:
        DESEQ2_DIR + "01-htseq/{sample}_counts.txt"
    threads: 2
    params:
        strand = rnaseq_snakefile_helper.strand_option_htseq(config["alignment_options"]["library_type"])
    log:
        DESEQ2_DIR + "01-htseq/.log/{sample}_htseq_per_sample.log"
    shell:
        "module purge && module load watermelon_dependencies && "
        "export MKL_NUM_THREADS={threads} && " #these exports throttle numpy processes
        "export NUMEXPR_NUM_THREADS={threads} && "
        "export OMP_NUM_THREADS={threads} && "
        "python -m HTSeq.scripts.count "
        "   -f bam "
        "   -s {params.strand} "
        "   -m intersection-nonempty "
        "   -q {input.bams} "
        "   {input.gtf} "
        "   > {output}.tmp "
        "   2>&1 | tee {log} && "
        "mv {output}.tmp {output} "

rule deseq2_htseq_merge:
    input:
        sample_checksum = CONFIG_CHECKSUMS_DIR + "config-samples.watermelon.md5",
        sample_count_files = expand(DESEQ2_DIR + "01-htseq/{sample}_counts.txt",
                                    sample=config[SAMPLES_KEY])
    output:
        counts_filename = DESEQ2_DIR + "01-htseq/htseq_merged.txt",
        stats_filename = DESEQ2_DIR + "01-htseq/htseq_stats.txt",
    params:
        input_dir = DESEQ2_DIR + "01-htseq/",
        counts_filename = DESEQ2_DIR + "01-htseq/htseq_merged.txt",
        stats_filename = DESEQ2_DIR + "01-htseq/htseq_stats.txt",
    log:
        DESEQ2_DIR + "01-htseq/.log/htseq_merge.log"
    shell:
       """(python {WATERMELON_SCRIPTS_DIR}/deseq2_htseq_merge.py \
       --htseq_dir={params.input_dir} \
       --suffix=_counts.txt \
       --counts_filename={output.counts_filename} \
       --stats_filename={output.stats_filename}
       ) 2>&1 | tee {log}"""


rule deseq2_metadata_contrasts:
    input:
        sample_checksum      = CONFIG_CHECKSUMS_DIR + "config-samples.watermelon.md5",
        comparison_checksum  = CONFIG_CHECKSUMS_DIR + "config-comparisons.watermelon.md5",
        phenotype_checksum   = CONFIG_CHECKSUMS_DIR + "config-phenotypes.watermelon.md5",
        main_factor_checksum = CONFIG_CHECKSUMS_DIR + "config-main_factors.watermelon.md5"
    output:
        sample_metadata = DESEQ2_DIR + "02-metadata_contrasts/sample_metadata.txt",
        contrasts = DESEQ2_DIR + "02-metadata_contrasts/contrasts.txt"
    params:
        phenos_with_replicates = phenotypeManager.phenotypes_with_replicates
    run:
        deseq2_helper.build_sample_metadata(config, params.phenos_with_replicates, output.sample_metadata)
        deseq2_helper.build_contrasts(config, params.phenos_with_replicates, output.contrasts)

rule deseq2_diffex:
    input:
        htseq_counts = DESEQ2_DIR + "01-htseq/htseq_merged.txt",
        sample_metadata = DESEQ2_DIR + "02-metadata_contrasts/sample_metadata.txt",
        contrasts = DESEQ2_DIR + "02-metadata_contrasts/contrasts.txt",
    output:
        dir = DESEQ2_DIR + "03-deseq2_diffex",
        files = expand(DESEQ2_DIR + "03-deseq2_diffex/gene_lists/{phenotype}/{comparison}.txt",
                       zip,
                       phenotype=REPLICATE_PHENOTYPE_NAMES,
                       comparison=REPLICATE_COMPARISON_GROUPS),
    threads: 8
    log:
        DESEQ2_DIR + "03-deseq2_diffex/.log/deseq2_DESeq2Diffex.log"
    params:
        fold_change = config["fold_change"],
        adjusted_pvalue = config["deseq2_adjustedPValue"],
    resources:
        memoryInGb = 16
    shell:
        "module purge && module load watermelon_dependencies && "
        "rm -rf {output.dir}/normalized_data && "
        "rm -rf {output.dir}/plots && "
        "rm -rf {output.dir}/gene_lists && "
        "rm -rf {output.dir}/.tmp/* && "
        "{WATERMELON_SCRIPTS_DIR}/deseq2_diffex.R "
        "    -c {input.htseq_counts} "
        "    -m {input.sample_metadata} "
        "    -f {input.contrasts} "
        "    -o {output.dir}/.tmp "
        "    --foldChange={params.fold_change} "
        "    --adjustedPValue={params.adjusted_pvalue} "
        "    --threads={threads} "
        "    --javaMemoryInGb={resources.memoryInGb} 2>&1 "
        "    --pandocMemoryInGb={resources.memoryInGb} 2>&1 | tee {log} && "
        "mv {output.dir}/.tmp/* {output.dir} && "
        "touch {output.dir} && "
        "rm -f Rplots.pdf " #Some part of R generates this empty (nuisance) plot

rule deseq2_annotation:
   input:
       diffex_file= DESEQ2_DIR + "03-deseq2_diffex/gene_lists/{phenotype}/{comparison}.txt",
       gene_info = "references/entrez_gene_info",
   output:
       DESEQ2_DIR + "04-annotation/{phenotype}/{comparison}.annot.txt",
   params:
       genome = config["genome"],
       output_dir = DESEQ2_DIR + "04-annotation/{comparison}",
   shell:
       "rm -rf {params.output_dir}.tmp && "
       "mkdir -p {params.output_dir}.tmp && "
       "python {WATERMELON_SCRIPTS_DIR}/deseq2_annotate.py "
       "    -i {input.gene_info} "
       "    -e {input.diffex_file} "
       "    -g {params.genome} "
       "    -o {output}.tmp && "
       "mv {output}.tmp {output} "

rule deseq2_run_info:
    input:
        glossary = WATERMELON_SCRIPTS_DIR + "deseq2_glossary.txt",
        sample_metadata = DESEQ2_DIR + "02-metadata_contrasts/sample_metadata.txt",
        contrasts = DESEQ2_DIR + "02-metadata_contrasts/contrasts.txt"
    output:
        run_info = DESEQ2_DIR + "05-run_info/run_info.txt",
        glossary = DESEQ2_DIR + "05-run_info/glossary.txt"
    run:
        command = ('module load watermelon_dependencies && '
                   'module list -t 2> {}').format(output['run_info'])
        subprocess.call(command, shell=True)
        with open(output['run_info'], 'a') as run_info_file:
            print('\n\nConfig\n', file=run_info_file)
            print(yaml.dump(config, default_flow_style=False),
                  file=run_info_file)
        copyfile(input['glossary'], output['glossary'])

rule deseq2_excel:
    input:
        gene = DESEQ2_DIR + "04-annotation/{phenotype}/{comparison}.annot.txt",
        glossary = DESEQ2_DIR + "05-run_info/glossary.txt",
        run_info = DESEQ2_DIR + "05-run_info/run_info.txt"
    output:
        annotated_file = DESEQ2_DIR + "06-excel/{phenotype}/{comparison}.xlsx",
    log:
        DESEQ2_DIR + "06-excel/.log/{phenotype}_diffex_excel.log"
    params:
        output_dir = DESEQ2_DIR + "06-excel/",
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py "
        " -g {input.gene}"
        " --glossary {input.glossary} "
        " --info_filepath {input.run_info} "
        " {output.annotated_file} "
        " 2>&1 | tee {log} "

rule deseq2_summary:
    input:
        input_files = expand(DESEQ2_DIR + "04-annotation/{phenotype}/{comparison}.annot.txt",
                             zip,
                             phenotype=REPLICATE_PHENOTYPE_NAMES,
                             comparison=REPLICATE_COMPARISON_GROUPS),
    output:
        summary_txt = DESEQ2_DIR + "07-summary/deseq2_summary.txt",
        summary_xlsx = DESEQ2_DIR + "07-summary/deseq2_summary.xlsx",
    log:
        DESEQ2_DIR + "07-summary/.log/deseq2_summary.log"
    params:
        output_dir = DESEQ2_DIR + "07-summary/",
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_summary.py "
        " --annotation_column gene_id"
        " --annotation_null . "
        " --diffex_call_column Call "
        " --diffex_call_pass YES "
        " --output_file {output.summary_txt} "
        " --output_xlsx {output.summary_xlsx} "
        " {input.input_files} "
        " 2>&1 | tee {log} && "
        "touch {params.output_dir}"

rule deliverables_deseq2:
    input:
        diffex_dir = DESEQ2_DIR + "03-deseq2_diffex",
        gene_lists = expand(DESEQ2_DIR + "06-excel/{phenotype_name}/{comparison}.xlsx",
                            zip,
                            phenotype_name=REPLICATE_PHENOTYPE_NAMES,
                            comparison=REPLICATE_COMPARISON_GROUPS),
        summary_txt = DESEQ2_DIR + "07-summary/deseq2_summary.txt",
        summary_xlsx = DESEQ2_DIR + "07-summary/deseq2_summary.xlsx",
    output:
        deliverables_dir = DELIVERABLES_DIR + "deseq2",
        summary_gene_lists = DELIVERABLES_DIR + "deseq2/gene_lists/deseq2_summary.txt",

    params:
        tmp_dir = DELIVERABLES_DIR + "deseq2.tmp",
        source_gene_list_dir = DESEQ2_DIR + "06-excel"
    shell:
        """rm -rf {output.deliverables_dir} {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        cp -r {input.diffex_dir}/counts {params.tmp_dir}
        cp -r {input.diffex_dir}/plots {params.tmp_dir}
        cp -r {params.source_gene_list_dir} {params.tmp_dir}/gene_lists
        cp {input.summary_txt} {params.tmp_dir}/gene_lists
        cp {input.summary_xlsx} {params.tmp_dir}/gene_lists
        mv {params.tmp_dir} {output.deliverables_dir} """
