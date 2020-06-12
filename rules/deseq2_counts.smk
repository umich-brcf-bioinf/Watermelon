'''This rule imports DESeq2 data from RSEM,
    outputs the raw counts and normalized counts as calculated by DESeq2
'''
rule deseq2_counts_from_tximport_rsem:
    input:
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results', sample=config[SAMPLES_KEY]),
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results', sample=config[SAMPLES_KEY])
    output:
        count_data_rda = DIFFEX_DIR + 'deseq2/counts/count_data.rda',
        raw = DIFFEX_DIR + 'deseq2/counts/deseq2_raw_counts.txt',
        norm = DIFFEX_DIR + 'deseq2/counts/deseq2_depth_normalized_counts.txt',
        rlog = DIFFEX_DIR + 'deseq2/counts/deseq2_rlog_normalized_counts.txt'
    log:
        DIFFEX_DIR + 'deseq2/counts/.log/DESeq2_counts.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    singularity: 'docker://umichbfxcore/wat_diffex:0.1.1'
    params:
        rsem_dir = ALIGNMENT_DIR + '04-rsem_star_align',
        snakemake_rdata = DIFFEX_DIR + 'deseq2/counts/.deseq2_counts_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_counts.R'
