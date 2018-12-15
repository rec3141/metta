#!/bin/bash

#input $1 is the directory containing the (corrected) reads to normalize before input into co-assembly
readdir=$1
threads=48
entropy=0.95
mindepth=0 #all

rm $readdir/all_merged_lowcov.fq.gz
rm $readdir/all_pairs_lowcov.fq.gz
rm $readdir/all_singles_lowcov.fq.gz

for prefix in `ls $readdir/trimmed_*.fastq.gz | xargs -i basename {} | cut -f2 -d'_' | sort -u `; do

#normal
    F1=($readdir/trimmed_$prefix*_R1*gz)
    R1=($readdir/trimmed_$prefix*_R2*gz)
    S1=($readdir/trimmed_$prefix*_R_*gz)
#    if [ ! -e $S1 ]; then touch $S1; fi

echo "$F1\n$R1\n$S1"

#first run bbduk to remove low-entropy sequences
#then concatenate all seqs
#then run bbnorm to normalize sequences to max 400x depth

#merge pairs to get longer seqs
~/apps/bbmap/bbmerge.sh threads=$threads in1=$F1 in2=$R1 out=$readdir/tmp.bbmerge_lowcov.fq outu1=$readdir/tmp.bbmerge1_lowcov.fq outu2=$readdir/tmp.bbmerge2_lowcov.fq
~/apps/bbmap/bbduk.sh threads=$threads in=$readdir/tmp.bbmerge_lowcov.fq int=f out=stdout.fq entropy=$entropy >> $readdir/all_merged_lowcov.fq

#entropy filter pairs
~/apps/bbmap/bbduk.sh threads=$threads in=$readdir/tmp.bbmerge1_lowcov.fq in2=$readdir/tmp.bbmerge2_lowcov.fq out=$readdir/tmp.bbduk_lowcov.fq entropy=$entropy
rm $readdir/tmp.bbmerge*_lowcov.fq

#bbnorm doesn't allow stdin; output interleaved
~/apps/bbmap/bbnorm.sh threads=$threads mindepth=$mindepth in=$readdir/tmp.bbduk_lowcov.fq out=stdout.fq int=t >> $readdir/all_pairs_lowcov.fq
rm $readdir/tmp.bbduk_lowcov.fq

#singles
~/apps/bbmap/bbduk.sh threads=$threads in=$S1 out=$readdir/tmp.bbduk_lowcov.fq entropy=$entropy
~/apps/bbmap/bbnorm.sh threads=$threads mindepth=$mindepth in=$readdir/tmp.bbduk_lowcov.fq out=stdout.fq >> $readdir/all_singles_lowcov.fq
rm $readdir/tmp.bbduk_lowcov.fq

done

gzip $readdir/all_pairs_lowcov.fq
gzip $readdir/all_singles_lowcov.fq
gzip $readdir/all_merged_lowcov.fq
