#!/bin/bash

poetry run python synth.py data/synthetic_data/PD3851a.vcf -c 200
poetry run python synth.py data/synthetic_data/PD3890a+.vcf -c 200
poetry run python synth.py data/synthetic_data/PD3904a+.vcf -c 200