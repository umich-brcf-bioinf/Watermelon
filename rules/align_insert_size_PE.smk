rule align_insert_size_PE:
    input:
        fasta_file = "references/bowtie2_index/genome.fa",
        bbmap_ref = "references/bbmap_ref",
        raw_fastq_R1 = ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_R1_PE.fastq.gz",
        raw_fastq_R2 = ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_R2_PE.fastq.gz",
    output:
        ALIGNMENT_DIR + "04-insert_size/{sample}_read_stats.txt"
    params:
        sample_id = "{sample}",
        downsample_readcount = config["insert_size"]["downsample_readcount"],
        pairlen = config["insert_size"]["pair_length"],
        temp_dir = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}",
        downsampled_fastq_r1 = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}/{sample}_R1.downsampled.fastq.gz",
        downsampled_fastq_r2 = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}/{sample}_R2.downsampled.fastq.gz",
        ihist_file = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}/{sample}_ihist.txt",
        lhist_file = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}/{sample}_lhist.txt",
        mapped_file = ALIGNMENT_DIR + "04-insert_size/.tmp/{sample}/{sample}.sam",
        output_dir = ALIGNMENT_DIR + "04-insert_size/",
    threads:
        8
    resources:
        memoryInGb=32
    log:
        ALIGNMENT_DIR + "04-insert_size/.log/{sample}_read_stats.log"
    shell:
        '''(module purge && module load watermelon_dependencies &&
        mkdir -p {params.temp_dir} &&
        JAVA_OPTS="-XX:ParallelGCThreads=2" &&
        reformat.sh t={threads} \
            sampleseed=42 \
            samplereadstarget={params.downsample_readcount} \
            in1={input.raw_fastq_R1} in2={input.raw_fastq_R2} \
            out={params.downsampled_fastq_r1} out2={params.downsampled_fastq_r2} &&
        bbmap.sh -Xmx{resources.memoryInGb}g t={threads} \
            semiperfectmode=t \
            requirecorrectstrand=f \
            pairlen={params.pairlen} \
            path={input.bbmap_ref} \
            in1={params.downsampled_fastq_r1} \
            in2={params.downsampled_fastq_r2} \
            ihist={params.ihist_file} \
            lhist={params.lhist_file} \
            out={params.mapped_file} &&
        module purge && module load python/3.6.1 &&
        python {WATERMELON_SCRIPTS_DIR}/read_stats_bbmap_single.py \
            --input_dir {params.temp_dir} \
            --sample_id {params.sample_id} \
            --output_dir {params.output_dir}
        ) 2>&1 | tee {log}'''
