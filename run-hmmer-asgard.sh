#!/bin/sh
# run hmmer on one node for each .hmm file in directory
# input $1 is amino acid fasta file .faa

# ./run-hmmer-asgard.sh scaffolds.1000bpp.faa > jobs_to_run.sh

seqfile=$1

if [ "$1" == "" ]; then 
	echo "provide input sequence .faa file"
	exit 1; 
fi

NODE=bio
THREADS=2

tmpfile=`mktemp -u -p. | cut -f3 -d'.'`

#for hmm in `cut -f1 FOAM-onto_rel1.tsv | cut -f1 -d' ' | sort -u | grep -v L1`; do 
#	if [ ! -e $hmm.hmm ]; then
#	./find-kos.sh $hmm
#	fi;

for hmmfile in `ls *.*.hmm | grep -v FOAM`; do 
	hmm=`basename $hmmfile .hmm`;
	basefile=`basename $seqfile .faa`;
	hmmd=`echo $hmm | cut -f1 -d'.'`
	if [ ! -d "$hmmd" ]; then 
	mkdir -p $hmmd
	mkdir $hmmd/txt
#	mkdir $hmmd/tsv
	mkdir $hmmd/pfam
#	mkdir $hmmd/aln
	fi;

if [ -e "$hmmd/tsv/$hmm.$basefile.tsv" ]; then continue; fi;
	echo $hmm

	#hmmsearch searches hmm profiles against a sequence database
	#hmmscan searches protein sequences against a hmm profile database
	#using hmmsearch because https://cryptogenomicon.org/2011/05/27/hmmscan-vs-hmmsearch-speed-the-numerology/

#quicker, run as jobarray
#        HMMCMD="~/apps/hmmer/binaries/hmmsearch -o $hmmd/txt/$hmm.$basefile.txt --tblout $hmmd/tsv/$hmm.$basefile.tsv --notextw --cpu $THREADS -Z 1000000 -E 0.1 $hmm.hmm $seqfile"
        HMMCMD="~/apps/hmmer/binaries/hmmsearch -o /dev/null --tblout $hmmd/tsv/$hmm.$basefile.tsv --notextw --cpu $THREADS -Z 1000000 -E 0.1 $hmm.hmm $seqfile"
#	echo $HMMCMD
	echo "$HMMCMD" >> jobs_to_run.$tmpfile

#slow
#	cd $hmm
#        HMMCMD="~/apps/hmmer/binaries/hmmsearch -o txt/$hmm.$basefile.txt --pfamtblout pfam/$hmm.$basefile.pfam --max -A aln/$hmm.$basefile.aln --tblout tsv/$hmm.$basefile.tsv --notextw --cpu $THREADS -Z 1000000 -E 0.1 ./../$hmm.hmm ./../$seqfile"
#       	sbatch --partition=$NODE --ntasks=$THREADS --tasks-per-node=$THREADS  <<-EOF
#!/bin/sh
#eval "$HMMCMD"
#EOF

#cd ..

done;



# submit jobarray to chinook
split -l 28 -d jobs_to_run.$tmpfile jobs_to_run.$tmpfile.
rm jobs_to_run.$tmpfile

for file in jobs_to_run.$tmpfile.*; do 
	NN=`grep -c '^' $file`
	sbatch --array=1-$NN ./run-jobarray.sh $file $NN; 
done
