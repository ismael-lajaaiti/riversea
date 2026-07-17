ARG R_VERSION=4.5.3

# Required by INLA.
FROM rocker/geospatial:${R_VERSION}

# Required by zen4R.
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    libxml2-dev \
    libsecret-1-dev \
    git \
  && rm -rf /var/lib/apt/lists/*


# Quarto.
ARG QUARTO_VERSION=1.9.36
RUN curl -LO https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb \
  && apt-get install -y ./quarto-${QUARTO_VERSION}-linux-amd64.deb \
  && rm quarto-${QUARTO_VERSION}-linux-amd64.deb \
  && quarto install tinytex --no-prompt

WORKDIR /project

ENV RENV_PATHS_LIBRARY=/renv/library
ENV RENV_PATHS_CACHE=/renv/cache

COPY renv.lock renv.lock
COPY renv/activate.R renv/activate.R
COPY .Rprofile .Rprofile

RUN R -e " \
  options(renv.config.repos.override = c( \
    CRAN     = 'https://packagemanager.posit.co/cran/latest', \
    INLA     = 'https://inla.r-inla-download.org/R/stable', \
    rOpenSci = 'https://ropensci.r-universe.dev' \
  )); \
  renv::restore()"

COPY . .

RUN R -e "renv::install('.', dependencies = FALSE)"

# Package cache/library are built as root; make them readable so an
# RStudio Server session (running as the non-root `rstudio` user) can
# still follow renv's symlinks into the cache.
RUN chmod -R a+rX /renv

CMD ["R"]
