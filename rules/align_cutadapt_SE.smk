rule align_cutadapt_SE:
    input:
        trimming_options_checksum = CONFIG_CHECKSUMS_DIR + "config-trimming_options.watermelon.md5",
        raw_fastq = ALIGNMENT_DIR + "01-raw_reads/{sample}_{read}_SE.fastq.gz",
    output:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_{read}_SE.fastq.gz",
    params:
        base_quality_5prime = config["trimming_options"]["base_quality_5prime"],
        base_quality_3prime = config["trimming_options"]["base_quality_3prime"],
        trim_length_5prime = config["trimming_options"]["trim_length_5prime"],
        trim_length_3prime = config["trimming_options"]["trim_length_3prime"],
        trimming_options = rnaseq_snakefile_helper.cutadapt_options(config["trimming_options"])
    log:
        ALIGNMENT_DIR + "02-cutadapt/.log/{sample}_cutadapt.log"
    run:
        if rnaseq_snakefile_helper.cutadapt_options(config["trimming_options"]):
            shell('''(module purge && module load watermelon_dependencies &&
                set -x &&
                cutadapt -q {params.base_quality_5prime},{params.base_quality_3prime} \
                    -u {params.trim_length_5prime} \
                    -u -{params.trim_length_3prime} \
                    --trim-n -m 20 \
                    -o {output}.tmp.gz \
                    {input.raw_fastq} &&
                mv {output}.tmp.gz {output}
                ) 2>&1 | tee {log} ''')
        else:
            shell('''(ln -sf ../../{input.raw_fastq} {output} &&
                echo No trimming done
                ) 2>&1 |tee {log} ''')
