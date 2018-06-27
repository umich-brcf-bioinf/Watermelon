rule tuxedo_run_info:
    input:
        TUXEDO_DIR + "07-split/last_split",
        glossary = WATERMELON_SCRIPTS_DIR + "tuxedo_glossary.txt"
    output:
        run_info=TUXEDO_DIR + "08-run_info/run_info.txt",
        glossary=TUXEDO_DIR + "08-run_info/glossary.txt"
    run:
        command = ('module load watermelon_dependencies && '
                   'module list -t 2> {}').format(output['run_info'])
        subprocess.call(command, shell=True)
        with open(output['run_info'], 'a') as run_info_file:
            print('\n\nConfig\n', file=run_info_file)
            print(yaml.dump(config, default_flow_style=False),
                  file=run_info_file)
        copyfile(input['glossary'], output['glossary'])
