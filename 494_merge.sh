#!/bin/bash
#this program merges sample-normalized reads, re-normalizes for spades assembly
#input $1 is normalized read directory
#input $2 is minimum kmer depth for bbnorm
threads=48

normdir=$1
mindepth=$2

#cat $normdir/norm_Plate*_pairs.fastq.gz > $normdir/all_pairs.fastq.gz
#cat $normdir/norm_Plate*_merged.fastq.gz > $normdir/all_merged.fastq.gz
#cat $normdir/norm_Plate*_singles.fastq.gz > $normdir/all_singles.fastq.gz

#takes about an hour per iteration
#for mindepth in 3 2 1 0; do

~/apps/bbmap/bbnorm.sh threads=$threads mindepth=$mindepth in=$normdir/all_pairs.fastq.gz out=$normdir/norm_$mindepth"_"pairs.fastq.gz int=t tbr=t lt=$mindepth
~/apps/bbmap/bbnorm.sh threads=$threads mindepth=$mindepth in=$normdir/all_merged.fastq.gz out=$normdir/norm_$mindepth"_"merged.fastq.gz int=t tbr=t lt=$mindepth
~/apps/bbmap/bbnorm.sh threads=$threads mindepth=$mindepth in=$normdir/all_singles.fastq.gz out=$normdir/norm_$mindepth"_"singles.fastq.gz int=t tbr=t lt=$mindepth

#done;


