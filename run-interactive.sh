#!/bin/bash
# starts an interactive shell on slurm
# input $1 tells which partition to run it on
# to specify a particular node add e.g. -w n145
srun -p $1 --nodes=1 --exclusive --pty /bin/bash
