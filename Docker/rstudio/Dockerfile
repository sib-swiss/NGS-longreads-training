FROM rocker/rstudio

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y libz-dev
RUN apt-get install -y libxml2-dev
RUN apt-get install -y libxt-dev

COPY install_packages.R ./

RUN Rscript ./install_packages.R
