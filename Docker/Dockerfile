FROM rocker/r-ver:4.2.3

MAINTAINER johannes.hausmann@tron-mainz.de

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

COPY ../DESCRIPTION /tmp/DESCRIPTION

RUN apt update && apt install --no-recommends -y \
    build-essentials \
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
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt clean && \
    apt autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

RUN wget -P $BIN "https://raw.githubusercontent.com/rocker-org/rocker-versioned2/R4.2.3/scripts/bin/install2.r" && chmod +x $BIN/install2.r

RUN install2.r --error --skipinstalled -n 8 \
    tidyverse \
    devtools \
    remotes \
    BiocManager \
    vroom

RUN export SPLICE2NEO_VERSION=$(grep 'Version:' /tmp/DESCRIPTION | awk '{print $2}' && \
	Rscript -e "library(remotes); remotes::install_git("https://github.com/TRON-Bioinformatics/splice2neo.git", ref = 'v${SPLICE2NEO_VERSION}')"


WORKDIR /

RUN apt autoremove -y

# Remove build dependencies from container
RUN apt remove -y \
    build-essentials \
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
    unixodbc-dev


ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]
