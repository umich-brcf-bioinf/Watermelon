rule ballgown_run_info:
    input:
        annot_files = expand(BALLGOWN_DIR + '02-annotate/{phenotype}/{comparison}_gene.annot.txt',
                             zip,
                             phenotype=ALL_PHENOTYPE_NAMES,
                             comparison=ALL_COMPARISON_GROUPS) +
                       expand(BALLGOWN_DIR + '02-annotate/{phenotype}/{comparison}_isoform.annot.txt',
                              zip,
                              phenotype=ALL_PHENOTYPE_NAMES,
                              comparison=ALL_COMPARISON_GROUPS),
        glossary = WATERMELON_SCRIPTS_DIR + 'ballgown_glossary.txt'
    output:
        run_info = BALLGOWN_DIR + '03-run_info/run_info.txt',
        glossary = BALLGOWN_DIR + '03-run_info/glossary.txt'
    run:
        command = ('module load watermelon_dependencies/{} && '
                   'module list -t 2> {}').format(WAT_VER, output['run_info'])
        subprocess.call(command, shell=True)
        with open(output['run_info'], 'a') as run_info_file:
            print('\n\nConfig\n', file=run_info_file)
            print(yaml.dump(config, default_flow_style=False),
                  file=run_info_file)
        copyfile(input['glossary'], output['glossary'])
