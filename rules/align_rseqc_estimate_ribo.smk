rule align_rseqc_estimate_ribo:
    input:
        genome_bam = ALIGNMENT_DIR + "04-rsem_star_align/{sample}.genome.bam"
    output:
        in_bam = temp(ALIGNMENT_DIR + "05-rseqc/{sample}_ribo_content.in.bam"),
        out_bam = temp(ALIGNMENT_DIR + "05-rseqc/{sample}_ribo_content.ex.bam"),
        junk_bam = temp(ALIGNMENT_DIR + "05-rseqc/{sample}_ribo_content.junk.bam"),
        counts = ALIGNMENT_DIR + "05-rseqc/{sample}_ribo_count.txt",
        percentages = ALIGNMENT_DIR + "05-rseqc/{sample}_ribo_pct.txt",
    params:
        project_name = config['report_info']['project_name'],
        bed_ref = config['references']['ribo_bed'],
        out_prefix = ALIGNMENT_DIR + "05-rseqc/{sample}_ribo_content"
    resources: time_min=600, mem_mb=4000
    container: ENV_INFO['rseqc']['image_str']
    log:
        JOB_LOG_DIR + "align_rseqc_estimate_ribo_{sample}.log"
    shell:
        '''(split_bam.py -i {input.genome_bam} -r {params.bed_ref} -o {params.out_prefix} > {output.counts} &&
{WATERMELON_SCRIPTS_DIR}split_bam_percentage.py -i {output.counts} -l "Ribosomal" "Non-ribosomal" "QCFail/Unmapped" -n {wildcards.sample} > {output.percentages}
) 2>&1 | tee {log} '''
