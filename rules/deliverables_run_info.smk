import glob
import yaml
import shutil

rule deliverables_run_info:
    input:
        custom_envs = glob.glob(WORKFLOW_BASEDIR + '/rules/envs/*/*.yaml'),
        biocontainer_envs = WORKFLOW_BASEDIR + "/envs/biocontainer_envs.yaml",
        samplesheet = config['samplesheet']
    output:
        versions = DELIVERABLES_DIR + "run_info/env_software_versions.yaml",
        samplesheet = DELIVERABLES_DIR + "run_info/" + os.path.basename(config['samplesheet'])
    params:
        project_name = config['report_info']['project_name']
    run:
        def transform_dict(env_dict):
            new_dict = {}
            swstrings = env_dict['dependencies']
            for string in swstrings:
                sw, ver = string.split("=")
                new_dict[sw] = ver
            return(new_dict)

        # Entry for watermelon pipeline version
        dep_dict = {'watermelon': WAT_VER}
        # Iterate over conda yamls
        for env_yaml in input.custom_envs:
            env_name = os.path.splitext(os.path.basename(env_yaml))[0]
            with open(env_yaml) as yfile:
                env_dict = yaml.load(yfile, Loader=yaml.SafeLoader)

            dep_dict[env_name] = transform_dict(env_dict)

        # Add the biocontainer envs
        with open(input.biocontainer_envs, "r") as efile:
            env_dict = yaml.load(efile, Loader=yaml.SafeLoader)
            # Transform to grab just the version strings
            for envir in env_dict:
                # Make nested to be congruent with the rest
                env = {envir: env_dict[envir]['ver_str']}
                env_dict[envir] = env
            dep_dict.update(env_dict)


        with open(output.versions, 'w') as ver_file:
            ver_file.write("#This file contains software version information in the format:\n" +
                "#environment:\n" +
                "#    software_package:\n" +
                "#        software_version\n")
            yaml.dump(dep_dict, ver_file, default_flow_style=False, indent=4)

        shutil.copy(input.samplesheet, output.samplesheet)
