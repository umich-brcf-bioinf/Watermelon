import yaml
import shutil

# VER_INFO is available through an import in the parent workflow
rule deliverables_run_info:
    input:
        samplesheet = config['samplesheet']
    output:
        versions = DELIVERABLES_DIR + "run_info/env_software_versions.yaml",
        samplesheet = DELIVERABLES_DIR + "run_info/" + os.path.basename(config['samplesheet'])
    params:
        project_name = config['report_info']['project_name']
    run:
        with open(output.versions, 'w') as ver_file:
            ver_file.write("#This file contains software version information in the format:\n" +
                "#environment:\n" +
                "#    software_package:\n" +
                "#        software_version\n")
            yaml.dump(VER_INFO, ver_file, default_flow_style=False, indent=4)

        shutil.copy(input.samplesheet, output.samplesheet)
