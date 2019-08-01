rule align_cutadapt_PE:
    input:
        raw_fastq_R1 = ALIGNMENT_DIR + "02-gz_reads/{sample}_R1_PE.fastq.gz",
        raw_fastq_R2 = ALIGNMENT_DIR + "02-gz_reads/{sample}_R2_PE.fastq.gz",
    output:
        R1 = ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_R1_PE.fastq.gz",
        R2 = ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_R2_PE.fastq.gz",
    params:
        base_quality_5prime = config["trimming_options"]["base_quality_5prime"],
        base_quality_3prime = config["trimming_options"]["base_quality_3prime"],
        trim_length_5prime = config["trimming_options"]["trim_length_5prime"],
        trim_length_3prime = config["trimming_options"]["trim_length_3prime"]
    log:
        ALIGNMENT_DIR + "02-cutadapt/.log/{sample}_cutadapt.log"
    conda:
        'envs/cutadapt.yaml'
    shell:
        '''(cutadapt -q {params.base_quality_5prime},{params.base_quality_3prime} \
                -u {params.trim_length_5prime} \
                -u -{params.trim_length_3prime} \
                --trim-n -m 20 \
                -o {output.R1}.tmp.gz \
                -p {output.R2}.tmp.gz \
                {input.raw_fastq_R1} {input.raw_fastq_R2}
            mv {output.R1}.tmp.gz {output.R1}
            mv {output.R2}.tmp.gz {output.R2}
            ) 2>&1 | tee {log} '''
