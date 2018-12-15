#/bin/bash
# this program makes the bbnorm jobarray

tmpfile=`mktemp -u -p. | cut -f3 -d'.'`

split -l 24 -d bbnorm-jobarray.sh bbnorm-jobarray.$tmpfile.
rm bbnorm-jobarray.sh

for file in bbnorm-jobarray.$tmpfile.*; do
        NN=`grep -c '^' $file`
        sbatch --array=1-$NN bin/492_run-bbnorm-jobarray.sh $file $NN;
done

