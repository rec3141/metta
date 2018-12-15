#/bin/bash
# this program makes the spades jobarray

tmpfile=`mktemp -u -p. | cut -f3 -d'.'`

split -l 24 -d spades-jobarray.sh spades-jobarray.$tmpfile.
rm spades-jobarray.sh

for file in spades-jobarray.$tmpfile.*; do
        NN=`grep -c '^' $file`
        sbatch --array=1-$NN bin/36_run-spades-jobarray.sh $file $NN;
done

