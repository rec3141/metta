#!/bin/bash
#input $1 is mega-assembly to map against
#input $2 is directory of reads to map
module load bio/SAMtools/1.5-pic-intel-2016b

assembly=$1 #to map against
f1=`dirname $assembly`
f1=`basename $f1`
dir=$2 #with fastq.gz to map
d1=`basename $dir`

outfile=bbmap_jobarray.sh
rm $outfile

outdir=mapped/$f1/$d1
mkdir -p $outdir

~/apps/bbmap/bbmap.sh ref=$assembly

	D1=`dirname $dir`
#	for prefix in `ls $dir/trimmed*fastq.gz | xargs -I{} basename {} | cut -f2 -d'_' | uniq`; do
	for prefix in `ls $dir/trimmed*fastq.gz | cut -f3 -d'_' | uniq`; do

		R1=`ls $dir/*$prefix*R1*fastq.gz`
		R2=`ls $dir/*$prefix*R2*fastq.gz`

		CMD="~/apps/bbmap/bbmap.sh in=$R1 in2=$R2 outm=$outdir/$prefix.mapped.bam bamscript=$outdir/$prefix.bs.sh; sh $outdir/$prefix.bs.sh; rm $outdir/$prefix.bs.sh; rm $outdir/$prefix.mapped.bam"
		
		echo $CMD >> $outfile
	done; 


