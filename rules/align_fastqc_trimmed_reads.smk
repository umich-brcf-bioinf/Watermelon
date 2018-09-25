rule align_fastqc_trimmed_reads:
    input:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_{read_endedness}.fastq.gz"
    output:
        ALIGNMENT_DIR + "03-fastqc_reads/{sample}_trimmed_{read_endedness}_fastqc.html"
    log:
        ALIGNMENT_DIR + "03-fastqc_reads/.log/{sample}_trimmed_{read_endedness}_fastqc.log"
    params:
        fastqc_dir = ALIGNMENT_DIR + "03-fastqc_reads"
    threads:
        # fastqc is not multithreaded, but Java spawns way too many processes,
        # so this keeps Snakemake from overruning the process limit.
        2
    shell:
        '''(module purge
        module load watermelon_dependencies/{WAT_VER}
        fastqc {input} -o {params.fastqc_dir}
        ) 2>&1 | tee {log}'''
