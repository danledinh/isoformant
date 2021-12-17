#!/bin/bash

ml pandoc/2.11.3.2 && \
ml texlive/2018-foss-2017a && \
source activate sphinx

pandoc \
-H float_adjustment.tex whitepaper.md \
-o whitepaper.pdf \
--filter=pandoc-fignos \
--filter=pandoc-eqnos \
--filter=pandoc-tablenos \
--citeproc

source deactivate
