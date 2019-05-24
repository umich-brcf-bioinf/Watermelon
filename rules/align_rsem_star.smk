from os.path import join

#_star_config = config['genome_reference']['star']
_star_config = {'genome_dir' : '/ccmb/BioinfCore/ActiveProjects/trsaari/sandbox/RSEM_ref'}

rule rsem_star_align:
    input:
        fastq_files = lambda wildcards: rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_{read_endedness}.fastq.gz",
                SAMPLE_READS,
                wildcards.sample),
        #This portion determines if a new reference must be created
        genomeParameters = join(\
                ALIGNMENT_DIR,
               '03-rsem_star_genome_generate',
               'genomeParameters.txt'
               ),
    output:
        join(ALIGNMENT_DIR , '04-rsem_star_align', '{sample}.genes.results'),
        join(ALIGNMENT_DIR , '04-rsem_star_align', '{sample}.isoforms.results'),
        join(ALIGNMENT_DIR , '04-rsem_star_align', '{sample}.genome.bam'),
        join(ALIGNMENT_DIR , '04-rsem_star_align', '{sample}.transcript.bam'),
    log:
        join(ALIGNMENT_DIR + '04-rsem_star_align', '.log', '{sample}.rsem_star_align.log')
    conda:
        '../envs/rsem_star.yaml'
    benchmark:
        ALIGNMENT_DIR + '/benchmarks/rsem_star_align.{sample}.benchmark.txt'
    threads: 12
    resources:
        mem_gb=30
    params:
        rsem_ref_base = join(\
            ALIGNMENT_DIR,
            '03-rsem_star_genome_generate',
            config['rsem_ref_prefix']
            ),
        outFileNamePrefix = join(ALIGNMENT_DIR, '04-rsem_star_align', '{sample}'),
        paired_end = lambda wildcards, input: '--paired-end' if rnaseq_snakefile_helper.detect_paired_end_bool(input.fastq_files) else '',
        strand_flag = rnaseq_snakefile_helper.strand_option_rsem(config)
    shell: '''(
STAR_PATH=$(dirname $(which STAR))
rsem-calculate-expression \
    {params.paired_end} \
    {params.strand_flag} \
    --star \
    --star-path $STAR_PATH \
    --star-gzipped-read-file \
    --star-output-genome-bam \
    -p {threads} \
    {input.fastq_files} \
    {params.rsem_ref_base} \
    {params.outFileNamePrefix}
mv {params.outFileNamePrefix}.STAR.genome.bam {params.outFileNamePrefix}.genome.bam
)2>&1 | tee {log}'''
