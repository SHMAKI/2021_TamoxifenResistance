#!/bin/bash
#PBS -l select=1:ncpus=10:mem=30gb
#PBS -N STAR_sr_all
#PBS -j oe
#PBS -m e
#PBS -J 1-141
#$ -S /bin/bash
#$ -cwd

# -J 1-x job number -> $PBS_ARRAY_INDEX
indexGenomeFile="/user1/tanpaku/okada/magi/ReferenceGenome/hg38_star/"
srFile="/user1/tanpaku/okada/magi/MCF-7_TamR/MCF7TamR_sr_trimmed"
outFile="/user1/tanpaku/okada/magi/MCF-7_TamR/TamR_sr_STAR"

SCRIPT_DIR=$(cd $(dirname $0); pwd)

file_all=(`find ${srFile} -name '*.fq.gz' -type f`)
echo ${#file_all[@]}

ArrNo=$(($PBS_ARRAY_INDEX - 1))
file=${file_all[${ArrNo}]}
echo "mapping from: ${file}"
echo "PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"
echo "arrayno.: ${ArrNo}"
file2=${file%_R1_trimmed.fq.gz*}
file3=${file2##*/}

# Mapping with STAR
 STAR --genomeDir ${indexGenomeFile} --runThreadN 10 --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix ${outFile}/${file3}_STAR_ \
       --outSAMstrandField intronMotif --readFilesCommand zcat \
       --readFilesIn ${file}
# Make bam index
 samtools index ${outFile}/${file3}_STAR_Aligned.sortedByCoord.out.bam
