FROM continuumio/miniconda3:23.5.2-0-alpine as build

RUN \
    conda install conda-libmamba-solver

RUN \
    echo -e > ~/.condarc '\
solver: libmamba\n\
channels:\n\
    - conda-forge\n\
    - bioconda\n\
    - defaults'

RUN \
    conda create -y -p /env/ \
        'blast' \
        'openjdk >=8'

RUN \
    mkdir -p /opt/cider/ && \
    wget -O /opt/cider/cider.jar 'https://github.com/hartwigmedical/hmftools/releases/download/cider-v1.0.3/cider_v1.0.3.jar'

RUN \
    conda clean -yaf

# Move Conda environment into standard BioContainers base image
FROM quay.io/bioconda/base-glibc-busybox-bash:2.1.0

COPY --from=build /env/ /env/
COPY --from=build /opt/cider/ /opt/cider/

ENV PATH="/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/env/lib/:${LD_LIBRARY_PATH}"
