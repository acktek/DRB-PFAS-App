# Use Rocker Shiny as base image
FROM rocker/shiny:4.4.1

# Install system dependencies required for R packages (especially sf, rgdal)
RUN apt-get update && apt-get install -y \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /srv/shiny-server

# Copy renv files first for better caching
COPY renv.lock renv.lock

# Install renv and restore packages
RUN R -e "install.packages('renv', repos='https://cran.rstudio.com/')" && \
    R -e "renv::restore()"

# Copy all app files
COPY app.R app.R
COPY *.csv ./
COPY *.shp *.shx *.dbf *.prj *.cpg *.sbn *.sbx *.xml ./
COPY www/ www/

# Expose Shiny port
EXPOSE 3838

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app.R', host='0.0.0.0', port=3838)"]
