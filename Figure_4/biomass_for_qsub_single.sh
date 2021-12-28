#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-20:1
#$ -m a
#$ -o /dev/null
#$ -e /home/smagi/qsub_logs

module load python/3.8.5
# arg1 : model module name
python python_worker_single.py $1 $SGE_TASK_ID
