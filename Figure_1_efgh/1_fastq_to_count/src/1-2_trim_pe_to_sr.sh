#!/bin/bash

SCRIPT_DIR=$(cd $(dirname $0); pwd)
mkdir TamR_PE_trimmed_sr36bp
cd /Volumes/HDPC-UT/MCF7TamR_pe
OUTPUT_DIR=../TamR_PE_trimmed_sr36bp

for file in `\find . -name  '*_R1.fastq.gz'`
do
#trim_galore
trim_galore --output_dir ${SCRIPT_DIR}/TamR_PE_trimmed_sr36bp --fastqc --trim1 --gzip --three_prime_clip_R1 65 ${file}

done
