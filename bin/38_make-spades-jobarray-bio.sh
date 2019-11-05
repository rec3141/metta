#/bin/bash
# this program makes the spades jobarray
# input $1 is the amount of memory in Gb

tmpfile=`mktemp -u -p. | cut -f3 -d'.'`

split -l 24 -d spades-jobarray.sh spades-jobarray.$tmpfile.
rm spades-jobarray.sh

for file in spades-jobarray.$tmpfile.*; do
        NN=`grep -c '^' $file`
        sbatch --array=1-$NN --mem $3"G" bin/37_run-spades-jobarray-bio.sh $file $NN;
done

