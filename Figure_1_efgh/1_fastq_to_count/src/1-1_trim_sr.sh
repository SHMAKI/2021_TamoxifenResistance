#!/bin/bash

SCRIPT_DIR=$(cd $(dirname $0); pwd)
mkdir TamR_sr_trimmed
cd /Volumes/HDPC-UT/MCF7TamR_sr
OUTPUT_DIR=../MCF7TamR_sr_trimmed

for file in `\find . -name  '*_R1.fastq.gz'`
do
#trim_galore
trim_galore --output_dir ${SCRIPT_DIR}/TamR_sr_trimmed --fastqc --trim1 --gzip ${file}

done
