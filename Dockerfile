FROM bioconductor/bioconductor_docker:RELEASE_3_19

RUN R -e "BiocManager::install(c('doParallel', 'illuminaHumanv2.db', 'illuminaHumanv4.db', 'limma', 'oligo', 'readr', 'GEOquery', 'illumniaio', 'tidyverse', 'broom'))"

RUN mkdir /app

#ADD GPL10558.tsv /
#ADD retrieveFromGEO.R /
#ADD normalizeBeadChip.R /

WORKDIR /app
