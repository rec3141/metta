#!/bin/bash

#input $1 is the directory containing the (corrected) reads to normalize before input into co-assembly
readdir=$1
threads=16
entropy=0.95
mindepth=3 #high coverage

rm $readdir/all_merged_highcov.fq.gz
rm $readdir/all_pairs_highcov.fq.gz
rm $readdir/all_singles_highcov.fq.gz

for prefix in `ls $readdir/trimmed_*.fastq.gz | xargs -i basename {} | cut -f2 -d'_' | sort -u `; do

#normal
    F1=($readdir/trimmed_$prefix"_"*_R1_*)
    R1=($readdir/trimmed_$prefix"_"*_R2_*)
    S1=($readdir/trimmed_$prefix"_"*_R_*)

echo "$F1\n$R1\n$S1"

#first run bbduk to remove low-entropy sequences
#then concatenate all seqs
#then run bbnorm to normalize sequences to max 100x depth

#merge pairs to get longer seqs
~/apps/bbmap/bbmerge.sh threads=$threads in1=$F1 in2=$R1 out=$readdir/tmphigh.bbmerge.fq outu1=$readdir/tmphigh.bbmerge1.fq outu2=$readdir/tmphigh.bbmerge2.fq
~/apps/bbmap/bbduk.sh threads=$threads in=$readdir/tmphigh.bbmerge.fq int=f out=stdout.fq entropy=$entropy >> $readdir/all_merged_highcov.fq

#entropy filter pairs
~/apps/bbmap/bbduk.sh threads=$threads in=$readdir/tmphigh.bbmerge1.fq in2=$readdir/tmphigh.bbmerge2.fq out=$readdir/tmphigh.bbduk.fq entropy=$entropy

rm $readdir/tmphigh.bbmerge*.fq

#bbnorm doesn't allow stdin; output interleaved
~/apps/bbmap/bbnorm.sh mindepth=$mindepth threads=$threads in=$readdir/tmphigh.bbduk.fq out=stdout.fq int=t >> $readdir/all_pairs_highcov.fq
rm $readdir/tmphigh.bbduk.fq

#singles
~/apps/bbmap/bbduk.sh threads=$threads in=$S1 out=$readdir/tmphigh.bbduk.fq entropy=$entropy
~/apps/bbmap/bbnorm.sh mindepth=$mindepth threads=$threads in=$readdir/tmphigh.bbduk.fq out=stdout.fq >> $readdir/all_singles_highcov.fq
rm $readdir/tmphigh.bbduk.fq

done

gzip $readdir/all_pairs_highcov.fq
gzip $readdir/all_singles_highcov.fq
gzip $readdir/all_merged_highcov.fq
