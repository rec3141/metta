#provide input fasta assembly
module load bio/SAMtools/1.5-pic-intel-2016b

assembly=$1 #to map against
directory=$2 #with fastq to map
filterfile=$3 #file with regex to select samples, e.g. "Plate1"
outfile=$4

echo "#!/bin/bash" > $outfile

for dir in $directory; do
	D1=`dirname $dir`
#	echo $dir;
#normal
	for file in `ls $dir/trimmed*fastq* | cut -f1-4 -d'_' | uniq | grep -f $filterfile`; do

#EBI
#	for file in `ls $dir/trimmed*fastq* | cut -f1-2 -d'_' | uniq | grep -f $filterfile`; do
#		echo $file

#normal
		R1=`ls $file*R1*fastq*`
		R2=`ls $file*R2*fastq*`

#EBI
#		R1=`ls $file*_1*fastq*`
#		R2=`ls $file*_2*fastq*`
		CMD="~/apps/bbmap/bbmap.sh in=$R1 in2=$R2 out=`basename $file`.mapped.bam ref=$assembly bamscript=bs.sh; sh bs.sh;"
		echo $CMD >> $outfile
#		echo "done mapping";
	done; 
done;


