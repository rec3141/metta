#!/bin/bash
NODE=bio
THREADS=7
MEM=214
module load bio/SAMtools/1.5-foss-2016b

sbatch --partition=$NODE --ntasks=$THREADS --tasks-per-node=$THREADS --mem=$MEM"G" map-list.sh    
#sbatch --partition=$NODE --ntasks=$THREADS --tasks-per-node=$THREADS --mem=$MEM"G" map-sort.sh    
