rule align_hisat2:
    input:
        alignment_options_checksum = CONFIG_CHECKSUMS_DIR + "config-alignment_options.watermelon.md5",
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        transcriptome_index = expand('references/hisat2_index/genome.{n}.ht2', n = range(1,9)),
        fastq_files = lambda wildcards: rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_{read_endedness}.fastq.gz",
                SAMPLE_READS,
                wildcards.sample),
    output:
        bam = ALIGNMENT_DIR + "04-hisat2/{sample}.bam",
        summary = ALIGNMENT_DIR + "04-hisat2/{sample}_align_summary.txt",
    benchmark:
        "benchmarks/{sample}.align_hisat2.benchmark.txt"
    conda:
        'envs/align_hisat2_stringtie.yaml'
    params:
        transcriptome_index = 'references/hisat2_index/genome',
        strand_flag = rnaseq_snakefile_helper.strand_option_hisat2(config),
        sample = '{sample}',
        fastq_flag = lambda wildcards, input: rnaseq_snakefile_helper.hisat_detect_paired_end(input.fastq_files)
    log:
        ALIGNMENT_DIR + "04-hisat2/.log/{sample}_hisat2.log"
    threads: 8
    shell:
        '''(hisat2 \
            -q \
            --max-intronlen 500000 \
            {params.strand_flag} \
            --threads {threads} \
            --new-summary \
            --summary-file {output.summary} \
            -x {params.transcriptome_index} \
            {params.fastq_flag} \
            | samtools view -b \
            | samtools sort -T {params.sample} -o {output.bam} -
        samtools index -b {output.bam}
        ) 2>&1 | tee {log}'''
