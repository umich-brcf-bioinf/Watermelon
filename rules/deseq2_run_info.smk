rule deseq2_run_info:
    input:
        glossary = WATERMELON_SCRIPTS_DIR + "deseq2_glossary.txt",
        sample_metadata = DESEQ2_DIR + "02-metadata_contrasts/sample_metadata.txt",
        contrasts = DESEQ2_DIR + "02-metadata_contrasts/contrasts.txt"
    output:
        run_info = DESEQ2_DIR + "05-run_info/run_info.txt",
        glossary = DESEQ2_DIR + "05-run_info/glossary.txt"
    run:
        command = ('module load watermelon_dependencies/{} && '
                   'module list -t 2> {}').format(_WAT_VER, output['run_info'])
        subprocess.call(command, shell=True)
        with open(output['run_info'], 'a') as run_info_file:
            print('\n\nConfig\n', file=run_info_file)
            print(yaml.dump(config, default_flow_style=False),
                  file=run_info_file)
        copyfile(input['glossary'], output['glossary'])
