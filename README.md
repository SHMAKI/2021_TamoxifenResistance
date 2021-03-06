# 2021_TamoxifenResistance

This repository contains the source codes for the analysis of the cellular experiment and sequencing data, and mathematical simulation used in the following paper:
* Magi, S., Ki, S., Ukai, M. et al. A combination approach of pseudotime analysis and mathematical modeling for understanding drug-resistant mechanisms. *Sci. Rep.* **11**, 18511 (2021). https://doi.org/10.1038/s41598-021-97887-z

## Figure 1-bcd, cell culture experiment
#### Requirements
    * R
      - tidyverse

## Figure_1-efgh
### /1_fastq_to_count
source codes for generating gene count matrix files from .fastq files deposited in [DDBJ Sequence Read Archive (DRA)](https://www.ddbj.nig.ac.jp/dra/index-e.html) ([DRA004349](https://ddbj.nig.ac.jp/resource/sra-submission/DRA004349))
#### Requirements
    * trim_galore
    * STAR (+ index file)
    * samtools
    * featureCounts (+ .gtf file)

### /2_DEG_analysis
source codes for differential expression analysis and enrichment analysis
#### Requirements
    * R 3.5.2
      - tidyverse
      - DESeq2 1.22.2
      - ggpubr
      - rtracklayer
      - amap
      - ggdendro
      - tidyheatmap

## Figure_2and3
Sequencing data were deposited into the [DRA](https://www.ddbj.nig.ac.jp/dra/index-e.html) ([DRA009126](https://ddbj.nig.ac.jp/resource/sra-submission/DRA009126)).  
Count data were deposited into [Genomic Expression Archive (GEA)](https://www.ddbj.nig.ac.jp/gea/index-e.html) ([E-GEAD-468](https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-468/)).
#### Requirements
    * R 3.5.2
      - tidyverse
      - Seurat 3.1.1
      - Monocle3 0.1.2
      - reticulate
      - ggsci
      - scales
      - ComplexHeatmap
      - tidyheatmap
    * gtf file (GRCh38.81)

## Figure_4
Mathematical simulation of cell growth were carried out by using [BioMASS](https://github.com/biomass-dev/biomass) with some modifications.

#### Requirements
    * python 3.7 or higher
      - numpy
      - pandas
      - matplotlib
      - seaborn
      - jupyter
      - importlib

source codes for parameter fitting:
  `qsub biomass_for_qsub_single.sh 2021_Magi_TamR`

.ipynb file for simulation and making figures:
  `Figure4.ipynb`

## Figure_5
#### Requirements
    * R
      - tidyverse
      - plyr
