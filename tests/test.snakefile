rule all:
    input:
        expand("{sample}.txt",
                sample=config["samples"]),

rule make_file:
    output:
        '{sample}.txt'
    params:
        content = lambda wildcards: config["samples"][wildcards.sample]
    shell:
        'echo {params.content} > {output}'