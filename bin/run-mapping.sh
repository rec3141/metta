#!/bin/bash
#generates jobs for jobarray.sh

rm bbmap_to_run.sh
#rm rep_to_run.sh

DIR=~/work/GEOTRACES/
i=1

for SUBDIR in $DIR/GEO*; do
  echo $SUBDIR
  S1=`basename $SUBDIR`
  for PREFIX in `ls $SUBDIR/trimmed*fastq.gz | grep -v KD | xargs -I{} basename {} | cut -f1-3 -d'_' | sort -u`; do

   for BIN in `ls bin* | cut -f1 -d'.' | cut -f1 -d'_' | sort -u`; do
#  for BIN in allbins; do
  if [ ! -e $BIN"_"$S1"_"$PREFIX"_mapped.fasta" ]; then  
   FILE1=`ls $SUBDIR/$PREFIX*R1*fastq.gz`; #trimmed_26-10_S10_L001_R1_001.fastq.gz
   FILE2=`ls $SUBDIR/$PREFIX*R2*fastq.gz`;
   BUILD=`echo $BIN | cut -b4-`
#   BUILD=999
   if [ $i -eq 0 ]; then
   BBCMD="~/apps/bbmap/bbmap.sh in=$FILE1 in2=$FILE2 outm=$BIN"_"$S1"_"$PREFIX"_mapped.fasta" ref=$BIN.contigs.fasta build=$BUILD -Xmx24g;"
   else
   BBCMD="~/apps/bbmap/bbmap.sh in=$FILE1 in2=$FILE2 outm=$BIN"_"$S1"_"$PREFIX"_mapped.fasta" build=$BUILD -Xmx24g;"
   fi
   echo "$BBCMD">>bbmap_to_run.sh
  fi
  #separate the forward and reverse reads
#  REPCMD="~/apps/bbmap/repair.sh in1=$BIN"_"$S1"_"$PREFIX"_mapped.fastq" out1=$BIN"_"$S1"_"$PREFIX"_mapped_R1.fastq" out2=$BIN"_"$S1"_"$PREFIX"_mapped_R2.fastq""
#  echo "$REPCMD">>rep_to_run.sh
  done;
i=1;
done;
   echo "complete $SUBDIR"
done;

