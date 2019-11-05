#!/bin/bash
#input $1 is directory of trimmed reads to find corrected reads for and link to reads/corrected/DIR

readdir=$1
r1=`basename $readdir`

mkdir -p reads/corrected/$r1

for PREFIX in `ls $readdir/trimmed_*.fastq.gz | xargs -i basename {} | cut -f2 -d'_' | sort -u `; do

if [ -d `pwd`/assembly/spades/spades_$PREFIX*/corrected ]; then	
	ln -s `pwd`/assembly/spades/spades_$PREFIX*/corrected/*.gz `pwd`/reads/corrected/$r1/
#	ls -l `pwd`/assembly/spades/spades_$PREFIX*/corrected/*.gz
fi;
done


