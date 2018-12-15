#/bin/bash
# this program runs prokka annotation on mega co-assembly
# then runs kaiju to get identifications

# DON'T RUN THIS IN PARALLEL; kaiju uses all the RAM
# input $1 is the co-assembly
# input $2 is the the metabat depth matrix

module load lang/Anaconda3/2.5.0
source activate `pwd`/prokka

module load bio/BLAST+/2.6.0-pic-intel-2016b-Python-2.7.12
export PATH=bin/:$PATH
export PERL5LIB=/home/recollins/work/nanobase/Arctic_metagenomes/prokka/lib/site_perl/5.26.2:/home/recollins/work/nanobase/Arctic_metagenomes/prokka/lib/5.26.2

prokka --help

assembly=$1
metabat=$2

kaijudir=annotation/kaiju
prokkadir=annotation/prokka
vizbindir=binning/vizbin
sketchdir=annotation/sketch
spadesdir=assembly/spades
kjdb=~/sfosdna/reference_dbs/kaiju

echo "start $assembly"

d1=`dirname $assembly`
f1=`basename $assembly .fasta`
echo $f1

# to make prokka happy
sed -E 's/(^>.{37}).*/\1/' $assembly > $d1/$f1.assembly.fasta

#run prokka
if [ ! -d "$prokkadir/prokka_$f1" ]; then 
	echo "running prokka"
	prokka --metagenome --outdir $prokkadir/prokka_$f1 --prefix $f1 --locustag $f1 --force --fast --rawproduct $d1/$f1.assembly.fasta
fi;

echo "done $assembly"




