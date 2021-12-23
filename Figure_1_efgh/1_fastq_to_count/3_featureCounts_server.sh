#!/bin/bash
#PBS -l select=1:ncpus=10:mem=30gb
#PBS -N featureCounts
#PBS -j oe
#PBS -m e

gtfFile="/user1/tanpaku/okada/magi/ReferenceGenome/Homo_sapience.GRCh38.96.gtf"

featureCounts -t exon -g gene_id
-a ${gtfFile}  \
-o counts_ensembl_from_sr.tsv /*.bam
