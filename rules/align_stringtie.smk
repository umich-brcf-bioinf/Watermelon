rule align_stringtie:
    input:
        reference_checksum = CONFIG_CHECKSUMS_DIR + 'config-references.watermelon.md5',
        bams = ALIGNMENT_DIR + '04-hisat2/{sample}.bam',
        gtf = 'references/gtf',
    output:
        gtf = ALIGNMENT_DIR + '06-stringtie/{sample}.gtf',
        ballgown = expand(ALIGNMENT_DIR + '06-stringtie/ballgown/{sample}/{bgown_prefixes}.ctab',
                          sample='{sample}',
                          bgown_prefixes=['e2t','e_data','i2t','i_data','t_data']),
        prep_input = ALIGNMENT_DIR + '06-stringtie/{sample}_prepDE_input.txt'
    log:
        ALIGNMENT_DIR + '06-stringtie/.log/{sample}.align_stringtie.log'
    benchmark:
        'benchmarks/{sample}.align_stringtie.benchmark.txt'
    conda:
        'envs/align_hisat2_stringtie.yaml'
    params:
        sample = '{sample}',
        strand_flag = rnaseq_snakefile_helper.strand_option_stringtie(config),
        ballgown_path = ALIGNMENT_DIR + '06-stringtie/ballgown/{sample}/'
    threads: 8
    shell:
        '''(
        echo {params.sample} {output.gtf} > {output.prep_input}

        stringtie {input.bams} \
            -p {threads} \
            -G {input.gtf} \
            -b {params.ballgown_path} \
            -e \
            {params.strand_flag} \
            -o {output.gtf}
        ) 2>&1 | tee {log}'''
