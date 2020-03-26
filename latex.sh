#!/bin/bash

# clean
rm -v main.aux
rm -v main.bbl
rm -v main.bcf
rm -v main.blg
rm -v main.log
rm -v main.out
rm -v main.pdf
rm -v main.run.xml
rm -v main.synctex.gz

# ensure bail out on error
set -eo pipefail

# setup environment
source binder/env.sh

# run pdflatex + biber + pdflatex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
biber main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
