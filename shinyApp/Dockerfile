# Use the official Shiny base image (R + Shiny Server preinstalled)
FROM rocker/shiny:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libhdf5-dev \
  gfortran \
  patch \
  libglpk-dev \
  && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN install2.r --error --skipinstalled --ncpus 4 \
  shiny \
  data.table \
  Matrix \
  DT \
  magrittr \
  ggplot2 \
  ggrepel \
  gridExtra \
  ggsci \
  readr \
  qs \
  pals \
  shinyhelper \
  purrr \
  ggdendro \
  hdf5r \
  BiocManager \
  dplyr \
  ggpubr \
  && R -q -e "BiocManager::install(c('tradeSeq', 'SingleCellExperiment', 'slingshot', 'SummarizedExperiment'))" \
  && rm -rf /tmp/downloaded_packages \
  && strip /usr/local/lib/R/site-library/*/libs/*.so

# Copy app files to the image
COPY . /srv/shiny-server/

# Copy the Shiny Server config
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Give shiny user ownership of the app files
RUN chown -R shiny:shiny /srv/shiny-server

# Expose the default Shiny port
EXPOSE 3838


# Run the app (Shiny Server will auto-serve it from /srv/shiny-server)
CMD ["/usr/bin/shiny-server"]
