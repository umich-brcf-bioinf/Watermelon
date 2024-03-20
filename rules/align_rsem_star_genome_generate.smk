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
        JOB_LOG_DIR + 'rsem_star_genome_generate.log'
    container: 'docker://umichbfxcore/rsem_star:0.1.1'
    benchmark:
        join(ALIGNMENT_DIR, 'benchmarks', 'rsem_star_genome_generate.benchmark.txt')
    resources: cpus_per_task=12, mem_mb=50000, runtime=360
    params:
        project_name = config['report_info']['project_name'],
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
    -p {resources.cpus_per_task} \
    {input.fasta} \
    {params.rsem_ref_base}
) 2>&1 | tee {log}'''
