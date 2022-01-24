def build_secondpass_cmd(wildcards):
    # Weird .fastq.gz.tmp.gz extension is due to cutadapt required automatic format extension
    firstpass_fq_R1 = R1 = ALIGNMENT_DIR + "02-cutadapt/{}_R1_trimmed.fastq.gz.tmp.gz".format(wildcards.sample)
    firstpass_fq_R2 = R2 = ALIGNMENT_DIR + "02-cutadapt/{}_R2_trimmed.fastq.gz.tmp.gz".format(wildcards.sample)
    secondpass_fq_R1 = R1 = ALIGNMENT_DIR + "02-cutadapt/{}_R1_trimmed.fastq.gz.tmp2.gz".format(wildcards.sample)
    secondpass_fq_R2 = R2 = ALIGNMENT_DIR + "02-cutadapt/{}_R2_trimmed.fastq.gz.tmp2.gz".format(wildcards.sample)
    if not config["trimming"]["secondpass_args"]: # Handle case when present but empty
        return ""
    else:
        cmd_fmt = "&&\ncutadapt --cores 0 {} -o {} -p {} {} {}"
        cmd = cmd_fmt.format(
            config["trimming"]["secondpass_args"],
            secondpass_fq_R1,
            secondpass_fq_R2,
            firstpass_fq_R1,
            firstpass_fq_R2
        )
        return cmd

rule align_cutadapt_PE:
    input:
        raw_fastq_R1 = ALIGNMENT_DIR + "02-gz_reads/{sample}_R1.fastq.gz",
        raw_fastq_R2 = ALIGNMENT_DIR + "02-gz_reads/{sample}_R2.fastq.gz",
    output:
        R1 = ALIGNMENT_DIR + "02-cutadapt/{sample}_R1_trimmed.fastq.gz",
        R2 = ALIGNMENT_DIR + "02-cutadapt/{sample}_R2_trimmed.fastq.gz",
    params:
        project_name = config["report_info"]["project_name"],
        cutadapt_args = config["trimming"]["cutadapt_args"],
        secondpass_cmd = build_secondpass_cmd if "secondpass_args" in config["trimming"] else ""
    log:
        JOB_LOG_DIR + "align_cutadapt_PE_{sample}.log"
    conda: 'envs/cutadapt/cutadapt.yaml'
    resources: time_min=300, cpus=8
    container: ENV_INFO['cutadapt']['image_str']
    shell:
        '''(cutadapt --cores {resources.cpus} {params.cutadapt_args} \
-o {output.R1}.tmp.gz \
-p {output.R2}.tmp.gz \
{input.raw_fastq_R1} {input.raw_fastq_R2} {params.secondpass_cmd}
if [[ -f {output.R1}.tmp2.gz ]]
    then mv {output.R1}.tmp2.gz {output.R1} &&
    mv {output.R2}.tmp2.gz {output.R2} &&
    rm {output.R1}.tmp.gz &&
    rm {output.R2}.tmp.gz &&
    echo "done w/ 2 of 2 passes"
else
    mv {output.R1}.tmp.gz {output.R1} &&
    mv {output.R2}.tmp.gz {output.R2} &&
    echo "done w/ 1 of 1 pass"
fi) 2>&1 | tee {log} '''
