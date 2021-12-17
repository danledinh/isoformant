#!/bin/bash

source activate isoformant

jupyter nbconvert \
--ClearMetadataPreprocessor.enabled=True \
--ClearOutput.enabled=True \
--to markdown \
README.ipynb

source deactivate
