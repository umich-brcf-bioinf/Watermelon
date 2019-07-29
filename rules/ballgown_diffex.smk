rule ballgown_diffex:
    input:
       expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results', sample=config[SAMPLES_KEY]),
       expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results', sample=config[SAMPLES_KEY])
    output:
        iso_results = DIFFEX_DIR + '{factor_name}/ballgown/ballgown_isoform_results.txt',
        rda = DIFFEX_DIR + '{factor_name}/ballgown/ballgown_data.rda',
    log:
        DIFFEX_DIR + '{factor_name}/ballgown/.log/ballgown_diffex.log'
    conda:
        'envs/diffex.yaml'
    params:
        rsem_dir = ALIGNMENT_DIR + '04-rsem_star_align/',
        snakemake_rdata = DIFFEX_DIR + '{factor_name}/ballgown/snakemake.rda' #TWS DEBUG
    script:
        '../scripts/ballgown_diffex.R'
