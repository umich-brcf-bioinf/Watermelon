rule align_cutadapt_SE:
    input:
        raw_fastq = ALIGNMENT_DIR + "02-gz_reads/{sample}_R{read}.fastq.gz",
    output:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_R{read}_trimmed.fastq.gz",
    params:
        cutadapt_args = config["trimming_options"]["cutadapt_args"]
    log:
        ALIGNMENT_DIR + "02-cutadapt/.log/{sample}_R{read}.align_cutadapt.log"
    conda: 'envs/cutadapt/cutadapt.yaml'
    singularity: 'docker://umichbfxcore/cutadapt'
    threads: 10
    shell:
        '''(cutadapt --cores {threads} \
                {params.cutadapt_args} \
                -o {output}.tmp.gz \
                {input.raw_fastq}
            mv {output}.tmp.gz {output}
            ) 2>&1 | tee {log} '''
