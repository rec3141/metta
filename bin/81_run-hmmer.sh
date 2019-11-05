#!/bin/sh
# run hmmer on one node for each .hmm file in directory
# input $1 is amino acid fasta file .faa

# ./run-hmmer.sh scaffolds.1000bpp.faa > hmmer_jobarray.sh

seqfile=$1
rm hmmer_jobarray.sh

THREADS=2
foamdir=annotation/foam

if [ "$1" == "" ]; then 
	echo "provide input sequence .faa file"
	exit 1; 
fi

#for hmm in `cut -f1 FOAM-onto_rel1.tsv | cut -f1 -d' ' | sort -u | grep -v L1`; do 
#	if [ ! -e $hmm.hmm ]; then
#	./find-kos.sh $hmm
#	fi;

for hmmfile in `ls $foamdir/*.*.hmm | grep -v FOAM`; do 
	echo $hmmfile
	hmm=`basename $hmmfile .hmm`;
	basefile=`basename $seqfile .faa`;
	hmmd=$foamdir/`echo $hmm | cut -f1 -d'.'`

	if [ ! -d "$hmmd" ]; then 
	mkdir -p $hmmd
	fi;
	mkdir -p $hmmd/txt
	mkdir -p $hmmd/tsv

if [ -e "$hmmd/tsv/$hmm.$basefile.tsv" ]; then continue; fi;
	echo $hmm

	#hmmsearch searches hmm profiles against a sequence database
	#hmmscan searches protein sequences against a hmm profile database
	#using hmmsearch because https://cryptogenomicon.org/2011/05/27/hmmscan-vs-hmmsearch-speed-the-numerology/

# quicker, run as jobarray
#        HMMCMD="~/apps/hmmer/binaries/hmmsearch -o $hmmd/txt/$hmm.$basefile.txt --tblout $hmmd/tsv/$hmm.$basefile.tsv --notextw --cpu $THREADS -Z 1000000 -E 0.1 $hmm.hmm $seqfile"
        HMMCMD="~/apps/hmmer/binaries/hmmsearch -o /dev/null --tblout $hmmd/tsv/$hmm.$basefile.tsv --notextw --cpu $THREADS -Z 1000000 -E 0.1 $foamdir/$hmm.hmm $seqfile"
	echo "$HMMCMD" >> hmmer_jobarray.sh

done;
