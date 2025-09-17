FROM docker/r-ver:4.2.3

ENV SRC /usr/local/src
ENV BIN /usr/local/bin


COPY DESCRIPTION /usr/local/src

RUN apt update && apt install -y \
    build-essential \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev \
    default-libmysqlclient-dev \
    libpq-dev \
    libsasl2-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    libxtst6 \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libbz2-dev \
    unixodbc-dev \
    tini \
    git \
    locales \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt clean && \
    apt autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

RUN install2.r --error --skipinstalled -n 8 \
    readr \
    tibble \
    dplyr \
    magrittr \
    argparse \
    stringr \
    tidyr \
    purrr \
    devtools \
    remotes \
    BiocManager \
    vroom \
    furrr


RUN Rscript -e 'library(remotes); remotes::install_version("Matrix", version = "1.6-5", repos = "http://cran.us.r-project.org", Ncpus=8)'
RUN Rscript -e 'library(remotes); remotes::install_version("MASS", version = "7.3-60", repos = "http://cran.us.r-project.org", Ncpus=8)'
RUN Rscript -e 'BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19"), Ncpus=8)'
RUN export SPLICE2NEO_VERSION=$(grep 'Version:' /usr/local/src/DESCRIPTION | awk '{print $2}') && \
Rscript -e "library(remotes); remotes::install_git('https://github.com/TRON-Bioinformatics/splice2neo.git', ref = 'v${SPLICE2NEO_VERSION}', Ncpus=8)"

RUN apt autoremove -y

# Remove build dependencies from container
RUN apt remove -y \
    build-essential \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev \
    default-libmysqlclient-dev \
    libpq-dev \
    libsasl2-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    libxtst6 \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    unixodbc-dev \
    libbz2-dev \
    git \
    locales

WORKDIR /
ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]
