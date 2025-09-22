# RNA-seq Analysis Pipeline
Marianna Mercadante
## Docker: R Environment for Data Analysis

To ensure a consistent and reproducible environment for data analysis, it is recommended to use a Docker container. Docker allows you to create a virtual environment that contains all the necessary software and dependencies. This helps avoid issues related to different versions of R and packages, ensuring that your analysis runs smoothly even years later.

### The Dockerfile

The Dockerfile is the blueprint for building the Docker image. It contains commands that install software, configure packages, and set up the environment. Each command creates a new layer in the Docker image.

1. Base Image

We start with a stable base: Ubuntu 20.04 LTS.

```{dockerfile}
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
```
`FROM ubuntu:20.04`: specifies the base image.

`ENV DEBIAN_FRONTEND=noninteractive`: prevents interactive prompts during package installation.

2. Install System Dependencies

R and many packages require system libraries and build tools:
```{dockerfile}
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common dirmngr gpg curl wget build-essential \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libharfbuzz-dev \
    libfribidi-dev make cmake gfortran libxt-dev liblapack-dev libblas-dev \
    zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev pandoc git sudo ca-certificates && \
    rm -rf /var/lib/apt/lists/*
```

* Installs essential tools for compilation (build-essential, cmake, make, gfortran).

* Installs libraries for graphics, fonts, and data processing (libpng, libtiff, libjpeg, etc.).

* Cleans up the apt cache to reduce image size.

3. Install R

We add the CRAN repository and install the latest version of R:
```{dockerfile}
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/cran.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" \
       | tee /etc/apt/sources.list.d/cran-r.list

RUN apt-get update && apt-get install -y --no-install-recommends r-base r-base-dev && \
    rm -rf /var/lib/apt/lists/*
```

* Adds the CRAN GPG key and repository.

* Installs R and the R development package (r-base-dev).

4. Install R Packages

We install commonly used CRAN packages and Bioconductor packages including ComplexHeatmap:
```{dockerfile}
RUN R -e "install.packages(c('ggplot2','VennDiagram','pheatmap','circlize','RColorBrewer','cluster','mclust','grid','reshape2','dplyr','data.table','BiocManager'), repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c('DESeq2','sva','GenomicRanges','rtracklayer','ComplexHeatmap'), ask=FALSE, update=FALSE)"
```
5. Add Non-Root User

To avoid permission issues inside the container, we create a non-root user:

```{dockerfile}
ARG USERNAME=containeruser
RUN useradd -u 1000 -m -s /bin/bash $USERNAME
USER $USERNAME
WORKDIR /home/$USERNAME
```
6. Copy Your Script and Set Entry Point

Finally, we copy the R script and set it to run automatically when the container starts:

```{dockerfile}
COPY dockerscript.R /home/containeruser/dockerscript.R
ENTRYPOINT ["Rscript", "/home/containeruser/dockerscript.R"]
```
### Building the Docker Image

To build the Docker image, place the Dockerfile in your working directory and run:

```{dockerfile}
docker build -t mio_progetto_r:latest .
```

`-t mio_progetto_r`:latest gives the image a name and tag.

`.` specifies that the Dockerfile is in the current directory.

It may take a few minutes as it downloads R and all required packages. Docker caches layers, so rebuilding later will be faster if the Dockerfile hasn't changed.

### Running a Docker Container

Once the image is built, you can run a container with your data and output directories mounted:
```{dockerfile}
docker run --rm -v "${PWD}/data:/data" -v "${PWD}/outputs:/results" mio_progetto_r:latest
```

`--rm`: automatically removes the container when stopped.

`-v "${PWD}/data:/data"`: mounts your local data folder into the container at /data.

`-v "${PWD}/outputs:/results"`: mounts your local outputs folder into the container at /results.

The R script dockerscript.R will automatically run inside the container and can read/write files from these mounted directories.
So now you have a structure like this one:
```{bash}
Local machine (your computer)              Docker container
─────────────────────────────             ─────────────────────────────
C:/Users/Marianna/Desktop/progetto/       (container root, e.g. /)
├── data/        <───────────────┐       /data/       # volume mapped to local data/
├── outputs/     <───────────────┘       /results/    # volume mapped to local outputs/
├── Dockerfile
├── README.md
```
## Analysis Description
The bioinformatics analysis focuses on identifying and characterizing differentially expressed genes in two RNA-seq datasets (GSE244485 and GSE244486) that show a clear batch effect.

The key steps of the analysis are as follows:
* **Gene Annotation**:
  * The `GTF file` (`Homo_sapiens.GRCh38.114.gtf.gz`) is imported to create a `GenomicRanges` object.
  * This object contains the genomic coordinates and annotations for all genes, which are used in downstream analyses.
* **Data Preparation**:
  * The raw count data from both datasets are combined into a single matrix.
  * Sample information (condition, batch) is loaded.
* **Batch Effect Correction**:
  * A Principal Component Analysis (PCA) is performed on the raw data to visualize the batch effect.
  * The batch effect is corrected using two different strategies:
    1. Including the batch as a variable in the analysis model (with `DESeq2`)
    2. Using the ComBat-seq algorithm to directly remove the batch effect from the data
  * A new PCA is run on the corrected data to verify the effectiveness of the batch removal.
* **Differential Expression (DE) Analysis**:
  * DE analysis is conducted using `DESeq2` on the data processed by both correction strategies.
  * The results (in terms of DEGs) from the two methods are compared, and their overlap is visualized with Venn diagrams.
* **Clustering and Heatmaps**:
  * Heatmaps are generated to visualize the expression profiles of genes unique to each correction method.
  * A dendrogram and a heatmap of samples based on the top 500 most variable genes are created to show how samples cluster based on their conditions.
* **Genomic Proximity**:
  * An additional analysis is performed to identify pairs of differentially expressed genes, one coding and one non-coding, that are located within ±5kb of each other.
  * Violin plots are generated to compare the expression levels of the identified coding and non-coding genes.

## Results and Outputs
The script generates several files and plots in the outputs directory:
 * `venn_infected_vs_mock.png`: Comparison of DE genes found by the two methods for the "infected" condition.
 * `venn_bystander_vs_mock.png`: Comparison of DE genes for the "bystander" condition.
 * `heatmap_*.png`: Heatmaps showing gene expression profiles.
 * `ARI_comparison_500genes_samples.csv`: A table with Adjusted Rand Index (ARI) values for various clustering metrics.
 * `ARI_boxplot_*.png`: Boxplots visualizing the distribution of ARI values.
 * `DE_coding_noncoding_pairs_investigation.csv`: A table listing pairs of differentially expressed coding and non-coding genes that are in close proximity.
 * `violinplot_coding_vs_noncoding_*.png`: Violin plots showing the expression levels of coding vs. non-coding genes.












