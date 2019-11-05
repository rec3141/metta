#!/bin/bash
# this program produces a jobarray for bbmapping reads to mega-assembly
#input $1 is the directory of corrected reads to map
#input $2 is the assembly to map to

rm bbmap_jobarray.sh

readdir=$1
assembly=$2

i=1

S1=`basename $readdir`
F1=`basename $assembly .fasta`

mkdir -p
#trimmed_Plate4-99-9G_S456_L001_R1_001.fastq.00.0_0.cor.fastq.gz
#trimmed_Plate4-99-9G_S456_L001_R2_001.fastq.00.0_0.cor.fastq.gz
#trimmed_Plate4-99-9G_S456_L001_R_unpaired.00.0_0.cor.fastq.gz

for prefix in `ls $dir/trimmed*fastq.gz | xargs -I{} basename {} | cut -f2 -d'_' | sort -u`; do

  if [ ! -e $F1"_"$S1"_"$prefix"_mapped.fasta" ]; then  
   FILE1=`ls $dir/trimmed*$prefix*R1*fastq.gz`; #trimmed_26-10_S10_L001_R1_001.fastq.gz
   FILE2=`ls $dir/trimmed*$prefix*R2*fastq.gz`;
   BUILD=1
   if [ $i -eq 0 ]; then
   BBCMD="~/apps/bbmap/bbmap.sh in=$FILE1 in2=$FILE2 outm=$F1"_"$S1"_"$prefix"_mapped.fasta" ref=$assembly build=$BUILD -Xmx24g;"
   else
   BBCMD="~/apps/bbmap/bbmap.sh in=$FILE1 in2=$FILE2 outm=$F1"_"$S1"_"$prefix"_mapped.fasta" build=$BUILD -Xmx24g;"
   fi
   echo "$BBCMD">>bbmap_jobarray.sh
  fi

i=1;
done;

echo "complete $dir"

