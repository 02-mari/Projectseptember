# Base image
FROM ubuntu:20.04

# Evita prompt interattivi durante l'installazione
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies di sistema
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common dirmngr gpg curl wget build-essential \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libharfbuzz-dev \
    libfribidi-dev make cmake gfortran libxt-dev liblapack-dev libblas-dev \
    zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev pandoc git sudo ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Aggiungi repository CRAN
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/cran.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" \
       | tee /etc/apt/sources.list.d/cran-r.list

# Installa R (ultima versione disponibile)
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base r-base-dev && \
    rm -rf /var/lib/apt/lists/*

# Installa pacchetti CRAN
RUN R -e "install.packages(c('ggplot2','VennDiagram','pheatmap','circlize','RColorBrewer','cluster','mclust','grid','reshape2','dplyr','data.table','BiocManager'), repos='https://cloud.r-project.org')"

# Installa pacchetti Bioconductor inclusi ComplexHeatmap
RUN R -e "BiocManager::install(c('DESeq2','sva','GenomicRanges','rtracklayer','ComplexHeatmap'), ask=FALSE, update=FALSE)"

# Aggiungi utente non-root
ARG USERNAME=containeruser
RUN useradd -u 1000 -m -s /bin/bash $USERNAME
USER $USERNAME
WORKDIR /home/$USERNAME

# Copia lo script nella home dell'utente
COPY dockerscript.R /home/containeruser/dockerscript.R

# Comando fisso: esegue lo script automaticamente
ENTRYPOINT ["Rscript", "/home/containeruser/dockerscript.R"]
