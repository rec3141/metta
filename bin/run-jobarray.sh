#!/bin/bash
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -p t1small

#run as sbatch --array=1-$2 ./run-jobarray.sh $1 $2
# input $1 is the list of jobs to run
# input $2 is the size of the array

# ${SLURM_ARRAY_TASK_ID} is a number

#this runs the ${SLURM_ARRAY_TASK_ID}'th line of $1
sed -n ${SLURM_ARRAY_TASK_ID}p $1 | bash
