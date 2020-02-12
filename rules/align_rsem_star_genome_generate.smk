from os.path import join

rule align_rsem_star_genome_generate:
    input:
        fasta = config['references']['fasta'],
        gtf = config['references']['gtf']
    output:
        parametersFiles = join(\
                ALIGNMENT_DIR,
                '04-rsem_star_genome_generate',
                'genomeParameters.txt'
                ),
    log:
        join(ALIGNMENT_DIR, '04-rsem_star_genome_generate', '.log', 'rsem_star_genome_generate.log')
    conda: 'envs/rsem_star/rsem_star.yaml'
    singularity: 'docker://umichbfxcore/rsem_star'
    benchmark:
        join(ALIGNMENT_DIR, 'benchmarks', 'rsem_star_genome_generate.benchmark.txt')
    threads: 12
    params:
        sjdbOverhang = int(config['alignment_options']['read_length']) - 1,
        rsem_ref_base = join(\
            ALIGNMENT_DIR,
            '04-rsem_star_genome_generate',
            config['genome']
            ),
    shell: '''
(
STAR_PATH=$(dirname $(which STAR))
rsem-prepare-reference \
    --gtf {input.gtf} \
    --star \
    --star-path $STAR_PATH \
    --star-sjdboverhang {params.sjdbOverhang} \
    -p {threads} \
    {input.fasta} \
    {params.rsem_ref_base}
) 2>&1 | tee {log}'''
