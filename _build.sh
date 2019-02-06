#!/bin/sh

set -ev

# Installs system libraries needed for R packages
sudo apt-get update && apt-get -y install "^libquantlib*"

Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::epub_book')"

