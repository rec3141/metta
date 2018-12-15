#!/bin/bash
#trim reads using bbduk
#input $1 is directory of raw reads to trim

DIR=$1
D1=`basename $DIR`
echo $D1

for PREFIX in `ls $DIR/*.fastq.gz | xargs -i basename {} | cut -f1 -d'_' | sort -u`; do
 echo $PREFIX

# normal
 F1=($DIR/$PREFIX"_"*R1*)
 R1=($DIR/$PREFIX"_"*R2*)

# for EBI
# F1=($DIR/$PREFIX"_"*1*)
# R1=($DIR/$PREFIX"_"*2*)

 F2=`basename $F1`
 R2=`basename $R1`

mkdir -p reads/trimmed/$D1

if [ ! -e reads/trimmed/$D1/trimmed_$R2 ]; then
 ~/apps/bbmap/bbduk.sh -Xmx1g in1=$F1 in2=$R1 out1=reads/trimmed/$D1/trimmed_$F2 out2=reads/trimmed/$D1/trimmed_$R2 ref=adapters ref=phix ktrim=r k=23 mink=11 hdist=1 tpe tbo #ow
fi;

done;
