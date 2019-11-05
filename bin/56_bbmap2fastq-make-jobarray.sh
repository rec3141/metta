#!/bin/bash
#input $1 is mega-assembly to map against
#input $2 is directory of reads to map
module load bio/SAMtools/1.5-pic-intel-2016b

assembly=$1 #to map against
f1=`basename $assembly .fasta`
dir=$2 #with fastq.gz to map
d1=`basename $dir`

outfile=bbmap_jobarray.sh
rm $outfile

outdir=mapped/$d1
mkdir -p $outdir

~/apps/bbmap/bbmap.sh in=$assembly ref=$assembly

echo "reference complete"

#	for prefix in `ls $dir/trimmed*fastq.gz | xargs -I{} basename {} | cut -f2 -d'_' | uniq`; do # super slow
	for prefix in `ls $dir/trimmed*fastq.gz | cut -f2 -d'_' | uniq`; do
#	for prefix in `ls $dir/*fastq.gz | cut -f1 -d'_' | uniq`; do

	echo $prefix

		R1=`ls $dir/*$prefix*R1*fastq.gz`
		R2=`ls $dir/*$prefix*R2*fastq.gz`

		CMD="~/apps/bbmap/bbmap.sh in="$R1" in2="$R2" outm="$outdir"/"$f1"_"$prefix"_R1.fastq.gz outm2="$outdir"/"$f1"_"$prefix"_R2.fastq.gz"
		
		echo $CMD >> $outfile
	done; 


