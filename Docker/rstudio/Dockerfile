FROM rocker/rstudio:4.2.0

# non interactive for installing dependencies
ENV DEBIAN_FRONTEND=noninteractive

# installing dependencies for DESeq2 and rmarkdown
RUN apt-get update
RUN apt-get install -y libz-dev
RUN apt-get install -y libxml2-dev
RUN apt-get install -y libxt-dev

# copy and run script to install packages
COPY install_packages.R ./

RUN Rscript ./install_packages.R
