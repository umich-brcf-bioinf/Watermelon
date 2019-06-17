#TWS - Maybe we can remove this rule entirely - doesn't seem like it has a real purpose
rule align_fastqc_align:
    input:
        ALIGNMENT_DIR + "04-rsem_star_align/{sample}.genome.bam"
    output:
        ALIGNMENT_DIR + "05-fastqc_align/{sample}_fastqc.html"
    params:
        fastqc_dir =  ALIGNMENT_DIR + "05-fastqc_align"
    log:
        ALIGNMENT_DIR + "05-fastqc_align/.log/{sample}_fastqc_align.log"
    conda:
        '../envs/fastqc.yaml'
    threads:
        # fastqc is not multithreaded, but Java spawns way too many processes,
        # so this keeps Snakemake from overruning the process limit.
        2
    shell:
        'fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log} '
