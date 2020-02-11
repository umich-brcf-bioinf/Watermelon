# fastqc environment

fastqc.yaml contains conda recipe, with explicitly defined software version

Dockerfile is generalized - env_name is passed in during docker build, and conda recipe is used to build the docker image
