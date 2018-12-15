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

kaijubin=~/apps/kaiju/bin
$kaijubin/kaiju

assembly=$1
metabat=$2

kaijudir=annotation/kaiju
prokkadir=annotation/prokka
kjdb=~/work/sfosdna/reference_dbs/kaiju

echo "start $assembly"

echo "not working? doesn't output taxonomy names in final files"

d1=`dirname $assembly`
f1=`basename $assembly .fasta`
echo $f1

# run-prokka.sh here

#run kaiju
if [ ! -e "$kaijudir/$f1.kaiju.tsv" ]; then
	echo "running kaiju"
	$kaijubin/kaiju -t $kjdb/nodes.dmp -f $kjdb/kaiju_db_nr_euk.fmi -p -i $prokkadir/prokka_$f1/$f1.faa -a greedy -e 5 -o $kaijudir/$f1.kaiju.tsv -z 24 -v
fi;
