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
        JOB_LOG_DIR + "align_cutadapt_PE_{sample}.log"
    conda: 'envs/cutadapt/cutadapt.yaml'
    singularity: 'docker://umichbfxcore/cutadapt'
    shell:
        '''(cutadapt {params.cutadapt_args} \
-o {output.R1}.tmp.gz \
-p {output.R2}.tmp.gz \
{input.raw_fastq_R1} {input.raw_fastq_R2}
mv {output.R1}.tmp.gz {output.R1} &&
mv {output.R2}.tmp.gz {output.R2}
) 2>&1 | tee {log} '''
