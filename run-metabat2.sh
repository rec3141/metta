#provide input fasta assembly
assembly=$1

for DIR in */; do
	D1=`dirname $DIR`
#	echo $DIR;
	for FILE in `ls $DIR/trimmed*fastq* | cut -f1-4 -d'_' | uniq | grep -f ice-samples2.txt`; do
#		echo $FILE
		R1=`ls $FILE*R1*fastq*`
		R2=`ls $FILE*R2*fastq*`
		CMD="~/apps/bbmap/bbmap.sh in=$R1 in2=$R2 out=`basename $FILE`.mapped.bam ref=$assembly bamscript=bs.sh; sh bs.sh;"
		echo $CMD 
#		echo "done mapping";
	done; 
done;

#/work/cryomics/apps/metabat/runMetaBat.sh $assembly *.mapped_sorted.bam;
#echo "done binning"; 

