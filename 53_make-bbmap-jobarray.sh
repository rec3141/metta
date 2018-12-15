#/bin/bash
# this program runs the bbmap jobarray
module load bio/SAMtools/1.5-pic-intel-2016b

tmpfile=`mktemp -u -p. | cut -f3 -d'.'`

split -l 24 -d bbmap_jobarray.sh bbmap_jobarray.$tmpfile.
rm bbmap_jobarray.sh

for file in bbmap_jobarray.$tmpfile.*; do
        NN=`grep -c '^' $file`
        sbatch --array=1-$NN bin/52_run-bbmap-jobarray.sh $file $NN;
done

