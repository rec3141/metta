#!/bin/bash
#input $1 is the directory containing the (corrected) reads to normalize before input into co-assembly
readdir=$1

rm bbnorm-jobarray.sh
for prefix in `ls $readdir/trimmed_*.fastq.gz | xargs -i basename {} | cut -f2 -d'_' | sort -u `; do

echo bin/493_run-bbmap-all-qc.sh $readdir $prefix >> bbnorm-jobarray.sh

done

bin/491_make-bbmap-jobarray.sh
