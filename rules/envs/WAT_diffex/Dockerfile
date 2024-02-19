FROM continuumio/miniconda3:23.10.0-1

ARG env_name

# env_name is supplied as --build-arg to docker, and is identical between yaml file basename and environment name specified within it
COPY ${env_name}.yaml /tmp/

RUN conda install mamba -n base -c conda-forge

RUN mamba env create -f /tmp/${env_name}.yaml && conda clean --all -y

ENV PATH /opt/conda/envs/${env_name}/bin:$PATH

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/${env_name}/lib

RUN apt-get update && \
    apt-get install -y texlive-latex-base texlive-latex-recommended texlive-fonts-recommended texlive-plain-generic

# Add RStudio to environment, to simplify the experience of interactive / manual interventions
RUN apt-get update && \
    apt-get install -y libnss3 libasound2 libegl1 && \
    wget -P /tmp/ https://s3.amazonaws.com/rstudio-ide-build/desktop/bionic/amd64/rstudio-2022.07.2-576-amd64.deb && \
    apt install -y /tmp/rstudio-2022.07.2-576-amd64.deb

