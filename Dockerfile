FROM bioconductor/bioconductor_docker:RELEASE_3_19

#RUN R -e "BiocManager::install('data.table')"
RUN R -e "BiocManager::install('doParallel')"
RUN R -e "BiocManager::install('illuminaHumanv2.db')"
RUN R -e "BiocManager::install('illuminaHumanv4.db')"
RUN R -e "BiocManager::install('limma')"
RUN R -e "BiocManager::install('oligo')"
RUN R -e "BiocManager::install('readr')"

ADD normalizeBeadChip.R /

WORKDIR /
