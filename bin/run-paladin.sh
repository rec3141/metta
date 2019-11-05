#!/bin/sh
# run paladin on forward and reverse reads
# input $1 is amino acid fasta file .faa
# input $2 is read directory
# beforehand set up paladin index
# paladin index -r3 plass-assembly2.faa

# run as ./run-paladin.sh plass-assembly2.faa ~/work/msl464/Arctic_metagenomes/corrected

module load bio/SAMtools/1.5-pic-intel-2016b

seqfile=$1
reads=$2

NODE=bio
THREADS=2

if [ "$1" == "" ]; then
	echo "provide input sequences"
	exit 1;
fi

tmpfile=`mktemp -u -p. | cut -f3 -d'.'`
rm $jobs_to_run.$tmpfile

for seqfile in `find $2 -name "*.cor.fastq.gz"`; do
	s1=`basename $seqfile .fastq.gz`
	PALCMD="if [ ! -e $s1.sorted.bam ]; then ~/apps/paladin/paladin align -t 4 -T 20 $1 $seqfile | samtools view -Sb - | samtools sort -o $s1.sorted.bam -; fi"
	echo "$PALCMD" >> jobs_to_run.$tmpfile
done;


split -l 7 -d jobs_to_run.$tmpfile jobs_to_run.$tmpfile.
rm jobs_to_run.$tmpfile

for file in jobs_to_run.$tmpfile.*; do
	NN=`grep -c '^' $file`
	sbatch --array=1-$NN ./run-jobarray.sh $file $NN;
done
