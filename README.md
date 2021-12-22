# 2021_TamoxifenResistance

This repository contains the source codes for the analysis of cellular experiment and sequencing data, and mathematical simulation used in the following paper:
* Magi S. et al., A combination approach of pseudotime analysis and mathematical modeling for understanding drug-resistant mechanisms. Sci Rep 11, 18511 (2021)
https://doi.org/10.1038/s41598-021-97887-z

## Figure 1-bcd, cell culture experiment
### Contents
  * .R files: codes for making figures
  * /data: experimental data
  * /Fig: generated figures
### Requirements
    * R 3.5.2
      - tidyverse 1.3.1

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
    - DESeq2 1.22.2
    - ggpubr 0.2.3
    - rtracklayer 1.42.2
    - amap 0.8-16
    - ggdendro 0.1-20
    - tidyheatmap 0.0.0.9000
    ***
    - org.Hs.eg.db 3.7.0
    - clusterProfiler 3.10.1
    ***
...
