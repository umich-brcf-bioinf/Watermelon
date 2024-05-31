import copy
import glob
import os
import yaml

_PIPE_VER = '2.8.4'

def _transform_conda_dict(env_dict):
    new_dict = {}
    swstrings = env_dict['dependencies']
    for string in swstrings:
        sw, ver = string.split("=")
        new_dict[sw] = ver
    return(new_dict)

# This snakefile gives access to VER_INFO (version info dict)
# and also to ENV_INFO (singularity info dict)
# Start with an entry for watermelon pipeline version
VER_INFO = {'watermelon': _PIPE_VER}
# Iterate over conda yamls
# Note WORKFLOW_BASEDIR is defined in parent snakefile
for env_yaml in glob.glob(WORKFLOW_BASEDIR + '/rules/envs/*/*.yaml'):
    env_name = os.path.splitext(os.path.basename(env_yaml))[0]
    with open(env_yaml) as yfile:
        env_dict = yaml.load(yfile, Loader=yaml.SafeLoader)

    VER_INFO[env_name] = _transform_conda_dict(env_dict)

ENV_INFO = {} # Context
# Add the biocontainer envs
with open(WORKFLOW_BASEDIR + "/envs/biocontainer_envs.yaml", "r") as efile:
    env_dict = yaml.load(efile, Loader=yaml.SafeLoader)
    ENV_INFO = copy.copy(env_dict) # Keep a copy of orig, make available to workflow
    # Transform to grab just the version strings
    for envir in env_dict:
        # Make nested to be congruent with the rest
        env = {envir: env_dict[envir]['ver_str']}
        env_dict[envir] = env
    VER_INFO.update(env_dict)
