ALL.append([DELIVERABLES_DIR + "combined_summary.txt", DELIVERABLES_DIR + "combined_summary.xlsx"])

def gene_summaries():
    summaries = [DELIVERABLES_DIR + "ballgown/gene_lists/ballgown_summary.txt"]
    if DESEQ2_ALL:
        summaries.append(DELIVERABLES_DIR + "deseq2/gene_lists/deseq2_summary.txt")
    return summaries

rule deliverables_combined_summary:
    input:
        lambda x: gene_summaries()
    output:
        txt = DELIVERABLES_DIR + "combined_summary.txt",
        xlsx = DELIVERABLES_DIR + "combined_summary.xlsx",
    params:
        base_name = DELIVERABLES_DIR + "combined_summary"
    shell:
        '''python {WATERMELON_SCRIPTS_DIR}/combine_summaries.py \
            --output_base {params.base_name} {input}
        '''
