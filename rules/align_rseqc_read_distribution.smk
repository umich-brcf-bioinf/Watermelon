rule align_rseqc_read_distribution:
    input:
        genome_bam = ALIGNMENT_DIR + "04-rsem_star_align/{sample}.genome.bam"
    output:
        read_dist = ALIGNMENT_DIR + "05-rseqc/{sample}_read_distribution.txt"
    params:
        project_name = config['report_info']['project_name'],
        bed_ref = config['references']['coverage_bed']
    resources: time_min=600, mem_mb=4000
    singularity: ENV_INFO['rseqc']['image_str']
    log:
        JOB_LOG_DIR + "align_rseqc_read_distribution_{sample}.log"
    shell:
        '''(read_distribution.py -i {input.genome_bam} -r {params.bed_ref} > {output.read_dist}.tmp
mv {output.read_dist}.tmp {output.read_dist}
) 2>&1 | tee {log} '''
