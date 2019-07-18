import glob
import yaml
import shutil

rule deliverables_run_info:
    input:
        envs = glob.glob(os.path.join(WORKFLOW_BASEDIR, 'rules', 'envs', '*.yaml')),
        conf = CONFIGFILE_PATH,
        samplesheet = config['sample_description_file']
    output:
        versions = DELIVERABLES_DIR + "run_info/env_software_versions.yaml",
        conf = DELIVERABLES_DIR + "run_info/" + os.path.basename(CONFIGFILE_PATH),
        samplesheet = DELIVERABLES_DIR + "run_info/" + os.path.basename(config['sample_description_file'])
    run:
        def transform_dict(env_dict):
            new_dict = {}
            swstrings = env_dict['dependencies']
            for string in swstrings:
                sw, ver = string.split("=")
                new_dict[sw] = ver
            return(new_dict)

        dep_dict = {}
        for env_yaml in input.envs:
            env_name = os.path.splitext(os.path.basename(env_yaml))[0]
            with open(env_yaml) as yfile:
                env_dict = yaml.load(yfile, Loader=yaml.SafeLoader)

            dep_dict[env_name] = transform_dict(env_dict)

        with open(output.versions, 'w') as ver_file:
            ver_file.write("#This file contains software version information in the format:\n" +
                "#environment:\n" +
                "#    software_package:\n" +
                "#        software_version\n")
            yaml.dump(dep_dict, ver_file, default_flow_style=False, indent=4)

        shutil.copy(input.conf, output.conf)
        shutil.copy(input.samplesheet, output.samplesheet)
