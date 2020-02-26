rule align_cutadapt_SE:
    input:
        raw_fastq = ALIGNMENT_DIR + "02-gz_reads/{sample}_R{read}.fastq.gz",
    output:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_R{read}_trimmed.fastq.gz",
    params:
        base_quality_5prime = config["trimming_options"]["base_quality_5prime"],
        base_quality_3prime = config["trimming_options"]["base_quality_3prime"],
        trim_length_5prime = config["trimming_options"]["trim_length_5prime"],
        trim_length_3prime = config["trimming_options"]["trim_length_3prime"]
    log:
        ALIGNMENT_DIR + "02-cutadapt/.log/{sample}_R{read}.align_cutadapt.log"
    conda: 'envs/cutadapt/cutadapt.yaml'
    singularity: 'docker://umichbfxcore/cutadapt'
    threads: 10
    shell:
        '''(cutadapt --cores {threads} \
                -q {params.base_quality_5prime},{params.base_quality_3prime} \
                -u {params.trim_length_5prime} \
                -u -{params.trim_length_3prime} \
                --trim-n -m 20 \
                -o {output}.tmp.gz \
                {input.raw_fastq}
            mv {output}.tmp.gz {output}
            ) 2>&1 | tee {log} '''
