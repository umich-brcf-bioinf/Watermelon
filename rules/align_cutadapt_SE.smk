def build_secondpass_cmd(wildcards):
    # Weird .fastq.gz.tmp.gz extension is due to cutadapt required automatic format extension
    firstpass_fq_R1 = R1 = ALIGNMENT_DIR + "02-cutadapt/{}_R1_trimmed.fastq.gz.tmp.gz".format(wildcards.sample)
    secondpass_fq_R1 = R1 = ALIGNMENT_DIR + "02-cutadapt/{}_R1_trimmed.fastq.gz.tmp2.gz".format(wildcards.sample)
    if not config["trimming"]["secondpass_args"]: # Handle case when present but empty
        return ""
    else:
        cmd_fmt = "&&\ncutadapt --cores 0 {} -o {} {}"
        cmd = cmd_fmt.format(
            config["trimming"]["secondpass_args"],
            secondpass_fq_R1,
            firstpass_fq_R1,
        )
        return cmd

rule align_cutadapt_SE:
    input:
        raw_fastq = ALIGNMENT_DIR + "02-gz_reads/{sample}_R{read}.fastq.gz",
    output:
        R1 = ALIGNMENT_DIR + "02-cutadapt/{sample}_R{read}_trimmed.fastq.gz",
    params:
        project_name = config['report_info']['project_name'],
        cutadapt_args = config["trimming"]["cutadapt_args"],
        secondpass_cmd = build_secondpass_cmd if "secondpass_args" in config["trimming"] else ""
    log:
        JOB_LOG_DIR + "align_cutadapt_SE_{sample}_R{read}.log"
    conda: 'envs/cutadapt/cutadapt.yaml'
    resources: time_min=300, cpus=8
    container: ENV_INFO['cutadapt']['image_str']
    shell:
        '''(cutadapt --cores {resources.cpus} \
{params.cutadapt_args} \
-o {output}.tmp.gz \
{input.raw_fastq} {params.secondpass_cmd}
if [[ -f {output.R1}.tmp2.gz ]]
    then mv {output.R1}.tmp2.gz {output.R1} &&
    rm {output.R1}.tmp.gz &&
    echo "done w/ 2 of 2 passes"
else
    mv {output.R1}.tmp.gz {output.R1} &&
    echo "done w/ 1 of 1 pass"
fi) 2>&1 | tee {log} '''
