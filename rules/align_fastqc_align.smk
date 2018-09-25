rule align_fastqc_align:
    input:
        ALIGNMENT_DIR + "04-hisat2/{sample}.bam"
    output:
        ALIGNMENT_DIR + "05-fastqc_align/{sample}_fastqc.html"
    params:
        fastqc_dir =  ALIGNMENT_DIR + "05-fastqc_align"
    log:
        ALIGNMENT_DIR + "05-fastqc_align/.log/{sample}_fastqc_align.log"
    threads:
        # fastqc is not multithreaded, but Java spawns way too many processes,
        # so this keeps Snakemake from overruning the process limit.
        2
    shell:
        '''(module purge
        module load watermelon_dependencies/{WAT_VER}
        fastqc {input} -o {params.fastqc_dir}
        ) 2>&1 | tee {log} '''
