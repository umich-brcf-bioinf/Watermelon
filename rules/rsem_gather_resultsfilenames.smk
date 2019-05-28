'''This rule gathers the sample names and corresponding results file names
into a file which can be used to import results into DESeq2 via tximport
'''

rule rsem_gather_resultsfilenames:
    input:
        gene_level = expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results', sample=config[SAMPLES_KEY]),
        transcript_level = expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results', sample=config[SAMPLES_KEY])
    output:
        gene_level = ALIGNMENT_DIR + '04-rsem_star_align/gene_results_files.txt',
        transcript_level = ALIGNMENT_DIR + '04-rsem_star_align/transcript_results_files.txt'
    log:
        ALIGNMENT_DIR + '04-rsem_star_align/.log/rsem_gather_resultsfilenames.log'
    shell:
        '''(
for file in {input.gene_level}
    do SNAME=$(basename $file ".genes.results")
    echo -e "$SNAME\t$file"
    done > {output.gene_level}
for file in {input.transcript_level}
    do SNAME=$(basename $file ".isoforms.results")
    echo -e "$SNAME\t$file"
    done > {output.transcript_level}
) 2>&1 | tee {log}
        '''
