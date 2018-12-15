#/bin/bash
# this program runs the hmmer jobarray

tmpfile=`mktemp -u -p. | cut -f3 -d'.'`

split -l 24 -d hmmer_jobarray.sh hmmer_jobarray.$tmpfile.
rm hmmer_jobarray.sh

for file in hmmer_jobarray.$tmpfile.*; do
        NN=`grep -c '^' $file`
        sbatch --array=1-$NN bin/82_run-hmmer-jobarray.sh $file $NN;
done

