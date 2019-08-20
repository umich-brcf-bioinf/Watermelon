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
    conda:
        'envs/rsem_star.yaml'
    benchmark:
        join(ALIGNMENT_DIR, 'benchmarks', 'rsem_star_genome_generate.benchmark.txt')
    threads: 12
    params:
        rsem_ref_base = join(\
            ALIGNMENT_DIR,
            '04-rsem_star_genome_generate',
            config['alignment_options']['rsem_ref_prefix']
            ),
    shell: '''
(
STAR_PATH=$(dirname $(which STAR))
rsem-prepare-reference \
    --gtf {input.gtf} \
    --star \
    --star-path $STAR_PATH \
    -p {threads} \
    {input.fasta} \
    {params.rsem_ref_base}
) 2>&1 | tee {log}'''