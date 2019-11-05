#!/bin/bash
#input $1 is directory to find trimmed reads
NODE=bio
TASKS=7
MEM=1400

dir=$1

spadesdir=assembly/spades
mkdir -p $spadesdir

rm spades-jobarray.sh

for prefix in `ls $dir/trimmed_*.fastq.gz | xargs -i basename {} | cut -f2 -d'_' | sort -u `; do

#normal
    F1=($dir/trimmed_$prefix"_"*R1*)
    R1=($dir/trimmed_$prefix"_"*R2*)

#EBI
#    F1=($dir/trimmed_$prefix"_"*1*)
#    R1=($dir/trimmed_$prefix"_"*2*)

    #echo "${F1[*]}"
    #echo "${R1[*]}"

#skip if already finished
    if [ -e $spadesdir/spades_$prefix*/scaffolds.fasta ]; then continue; fi;

##skip of already started (wait for cleanup to restart)
#    if [ -d $spadesdir/spades_$prefix*/K21 ]; then continue; fi;

#if already corrected, then continue
    if [ -d $spadesdir/spades_$prefix*/corrected ]; then
      outdir="spades_$prefix*"
      SPADESCMD="~/apps/SPAdes-3.13.0-Linux/bin/metaspades.py -o $spadesdir/$outdir -m $MEM -t $TASKS --restart-from last"
    else 
      outdir="spades_"$prefix"_"`date +'%F-%H-%M-%S'`
    SPADESCMD="~/apps/SPAdes-3.13.0-Linux/bin/metaspades.py -o $spadesdir/$outdir -m $MEM -t $TASKS" # --tmp-dir /tmp"

    #this is for combining reads from same sample on multiple runs
    J=$((${#F1[@]}-1))
    for I in `seq 0 $J`; do
        K=$((I+1))
        SPADESCMD=$SPADESCMD" --pe$K-1 ${F1[$I]} --pe$K-2 ${R1[$I]} "
    done

    fi;

	#delete temp folders if finished
	SPADESCMD=$SPADESCMD"; if [ -e $spadesdir/$outdir/scaffolds.fasta ]; then rm -rf $spadesdir/$outdir/K*; fi"
        echo -e $SPADESCMD >> spades-jobarray.sh

	continue #for testing without running on cluster or using job-array

        sbatch --partition=$NODE --nodes=1 --ntasks=$TASKS --tasks-per-node=$TASKS --mem=$MEM"G" <<-EOF
#!/bin/bash
eval "$SPADESCMD"
EOF
    
done


echo bin/38_make-spades-jobarray-bio.sh $MEM
