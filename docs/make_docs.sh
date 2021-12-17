#!/bin/bash

source activate sphinx

sphinx-apidoc -o ./source ../isoformant
make clean
make html

source deactivate