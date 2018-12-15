#!/bin/bash

module load lang/Python/2.7.12-pic-intel-2016b
~/apps/redundans/redundans.py -v -i mapped/cyanos/cyanos_R1.fastq.gz mapped/cyanos/cyanos_R2.fastq.gz mapped/cyanos/cyanos_R0.fastq.gz -f assembly/spades/spades_cyanos/scaffolds.fasta -o assembly/redundans/redundans_cyanos -t 14 --minLength 1


