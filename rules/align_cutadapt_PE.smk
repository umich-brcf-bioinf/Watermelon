rule align_cutadapt_PE:
    input:
        raw_fastq_R1 = ALIGNMENT_DIR + "02-gz_reads/{sample}_R1.fastq.gz",
        raw_fastq_R2 = ALIGNMENT_DIR + "02-gz_reads/{sample}_R2.fastq.gz",
    output:
        R1 = ALIGNMENT_DIR + "02-cutadapt/{sample}_R1_trimmed.fastq.gz",
        R2 = ALIGNMENT_DIR + "02-cutadapt/{sample}_R2_trimmed.fastq.gz",
    params:
        cutadapt_args = config["trimming_options"]["cutadapt_args"]
    log:
        ALIGNMENT_DIR + "02-cutadapt/.log/{sample}_cutadapt.log"
    conda: 'envs/cutadapt/cutadapt.yaml'
    singularity: 'docker://umichbfxcore/cutadapt'
    shell:
        '''(cutadapt {params.cutadapt_args} \
-o {output.R1} \
-p {output.R2} \
{input.raw_fastq_R1} {input.raw_fastq_R2}
) 2>&1 | tee {log} '''
