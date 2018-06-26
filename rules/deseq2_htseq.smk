rule deseq2_htseq:
    input:
        reference_checksum = CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        bams = ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_accepted_hits.bam",
        gtf = "references/gtf"
    output:
        DESEQ2_DIR + "01-htseq/{sample}_counts.txt"
    threads: 2
    params:
        strand = rnaseq_snakefile_helper.strand_option_htseq(config["alignment_options"]["library_type"])
    log:
        DESEQ2_DIR + "01-htseq/.log/{sample}_htseq_per_sample.log"
    shell:
        "module purge && module load watermelon_dependencies && "
        "export MKL_NUM_THREADS={threads} && " #these exports throttle numpy processes
        "export NUMEXPR_NUM_THREADS={threads} && "
        "export OMP_NUM_THREADS={threads} && "
        "python -m HTSeq.scripts.count "
        "   -f bam "
        "   -s {params.strand} "
        "   -m intersection-nonempty "
        "   -q {input.bams} "
        "   {input.gtf} "
        "   > {output}.tmp "
        "   2>&1 | tee {log} && "
        "mv {output}.tmp {output} "
