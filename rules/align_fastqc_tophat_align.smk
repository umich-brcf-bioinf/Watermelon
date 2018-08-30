rule align_fastqc_tophat_align:
    input:
        ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_accepted_hits.bam"
    output:
        ALIGNMENT_DIR + "05-fastqc_align/{sample}_accepted_hits_fastqc.html"
    params:
        fastqc_dir =  ALIGNMENT_DIR + "05-fastqc_align"
    log:
        ALIGNMENT_DIR + "05-fastqc_align/.log/{sample}_fastqc_tophat_align.log"
    shell:
        '''(module purge
        module load watermelon_dependencies/{WAT_VER}
        fastqc {input} -o {params.fastqc_dir}
        ) 2>&1 | tee {log} '''
