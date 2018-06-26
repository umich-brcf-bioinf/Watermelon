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
