ARG BASE_CONTAINER=jupyter/minimal-notebook:latest
FROM $BASE_CONTAINER

USER root

# update current packages list and install gcc from build-essential for genomeview
RUN sudo apt-get update && sudo apt-get -y install build-essential

# create conda environment
COPY --chown=${NB_UID}:${NB_GID} environment.yml /home/${NB_USER}/reqs/
RUN cd /home/${NB_USER}/reqs/ && \
    conda env update -n base -f environment.yml && \
    conda clean --all -f -y
    
# fix permissions so directory is user-writable
RUN fix-permissions "${CONDA_DIR}"
RUN fix-permissions "/home/${NB_USER}"
    
# adding additional pip installs to conda env
# COPY --chown=${NB_UID}:${NB_GID} pip_requirements.txt /home/${NB_USER}/reqs/
# RUN cd /home/${NB_USER}/reqs/ && \
# pip install -r pip_requirements.txt && \
# pip install -U git+https://github.com/nspies/genomeview.git@c1ac797642b86f8f411f75f5de68c4a9e6519078

# switch back to jovyan to avoid accidental container runs as root
USER ${NB_UID}