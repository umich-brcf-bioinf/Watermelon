rule align_create_transcriptome_index:
    input:
        alignment_options_checksum =  CONFIG_CHECKSUMS_DIR + "config-alignment_options.watermelon.md5",
        reference_checksum =  CONFIG_CHECKSUMS_DIR + "config-references.watermelon.md5",
        gtf = "references/gtf",
        bowtie2_index_dir = "references/bowtie2_index"
    output:
        ALIGNMENT_DIR + "04-tophat/transcriptome_index/transcriptome.fa"
    params:
        transcriptome_dir = "transcriptome_index",
        temp_dir =  ALIGNMENT_DIR + "04-tophat/.tmp",
        output_dir = ALIGNMENT_DIR + "04-tophat",
        strand = rnaseq_snakefile_helper.strand_option_tophat(config["alignment_options"]["library_type"])
    log:
        ALIGNMENT_DIR + "04-tophat/.log/create_transcriptome_index.log"
    shell:
        '''(module purge
        module load watermelon_dependencies/{WAT_VER}
        mkdir -p {params.temp_dir}
        rm -rf {params.temp_dir}/*
        tophat -G {input.gtf} \
            --library-type {params.strand} \
            --transcriptome-index={params.temp_dir}/transcriptome_index/transcriptome \
            {input.bowtie2_index_dir}/genome
        rm -rf {params.output_dir}/{params.transcriptome_dir}
        mv {params.temp_dir}/{params.transcriptome_dir} {params.output_dir}
        mv tophat_out {params.output_dir}/{params.transcriptome_dir}/
        touch {params.output_dir}/{params.transcriptome_dir}/*
        ) 2>&1 | tee {log}'''
