
rule align_rsem_star:
    input:
        fastq_files = lambda wildcards: expand(ALIGNMENT_DIR + "02-cutadapt/{basename}_trimmed.fastq.gz",
            basename=INPUT_MANAGER.sample_bnames_dict[wildcards.sample]
        ),
        #This portion determines if a new reference must be created
        genomeParameters = ALIGNMENT_DIR + '04-rsem_star_genome_generate/genomeParameters.txt'
    output:
        ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results',
        ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results',
        ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genome.bam',
        ALIGNMENT_DIR + '04-rsem_star_align/{sample}.transcript.bam',
        ALIGNMENT_DIR + '04-rsem_star_align/{sample}.stat/{sample}.cnt',
        ALIGNMENT_DIR + '04-rsem_star_align/{sample}.log'

    log:
        JOB_LOG_DIR + 'rsem_star_align_{sample}.log'
    conda: 'envs/rsem_star/rsem_star.yaml'
    singularity: 'docker://umichbfxcore/rsem_star:0.1.1'
    benchmark:
        ALIGNMENT_DIR + 'benchmarks/rsem_star_align.{sample}.benchmark.txt'
    resources: cpus=12, mem_mb=40000, time_min=720
    params:
        rsem_ref_base = ALIGNMENT_DIR + '04-rsem_star_genome_generate/' + config['genome'],
        outFileNamePrefix = ALIGNMENT_DIR + '04-rsem_star_align/{sample}',
        paired_end = lambda wildcards, input: '--paired-end' if rnaseq_snakefile_helper.detect_paired_end_bool(input.fastq_files) else ''
    shell: '''(
STAR_PATH=$(dirname $(which STAR))
rsem-calculate-expression \
    {params.paired_end} \
    --star \
    --star-path $STAR_PATH \
    --star-gzipped-read-file \
    --star-output-genome-bam \
    -p {resources.cpus} \
    {input.fastq_files} \
    {params.rsem_ref_base} \
    {params.outFileNamePrefix}
samtools sort -@ {resources.cpus} \
    -o {params.outFileNamePrefix}.genome.bam \
    {params.outFileNamePrefix}.STAR.genome.bam
rm {params.outFileNamePrefix}.STAR.genome.bam
)2>&1 | tee {log}'''
