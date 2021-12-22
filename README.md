# 2021_TamoxifenResistance

This repository contains the source codes for the analysis of cellular experiment and sequencing data, and mathematical simulation used in the following paper:
* Magi S. et al., A combination approach of pseudotime analysis and mathematical modeling for understanding drug-resistant mechanisms. Sci Rep 11, 18511 (2021)
https://doi.org/10.1038/s41598-021-97887-z

## Figure 1-bcd, cell culture experiment
### Requirement
  * R 3.5.2
    - tidyverse 1.3.1
### Contents
  * R files: codes for making figure
  * /data: experimental data
  * /Fig: generated figures

## Figure 1-efgh, bulk RNA-seq data
### 1_fastq_to_count
source codes for generating gene count matrix files from .fastq files ([DRA004349](https://ddbj.nig.ac.jp/resource/sra-submission/DRA004349))
#### Requirement
  * trim_galore
  * STAR
  * samtools
  * featureCounts

### 2_analysis
source codes for differential expression analysis and enrichment analysis
#### Requirement
  * R 3.5.2
    - tidyverse 1.3.1
...
