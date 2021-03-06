#!/bin/bash
#input $1 is directory to find trimmed reads
NODE=t2standard
THREADS=48
CPUS=24
MEM=120

DIR=$1

spadesdir=assembly/spades

for PREFIX in `ls $DIR/trimmed_*.fastq.gz | xargs -i basename {} | cut -f2 -d'_' | sort -u `; do

#normal
    F1=($DIR/trimmed_$PREFIX"_"*R1*)
    R1=($DIR/trimmed_$PREFIX"_"*R2*)

#EBI
#    F1=($DIR/trimmed_$PREFIX"_"*1*)
#    R1=($DIR/trimmed_$PREFIX"_"*2*)

    #echo "${F1[*]}"
    #echo "${R1[*]}"

    FINDIR=`basename $F1 | cut -f2 -d'_'` #isn't this just $PREFIX?
    OUTDIR=$FINDIR"_"`date +'%F-%H-%M-%S'`

    if [ -e assembly/spades/spades_$FINDIR*/scaffolds.fasta ]; then continue; fi;

    SPADESCMD="~/apps/SPAdes-3.13.0-Linux/bin/spades.py -o $spadesdir/spades_$OUTDIR -t $THREADS -m $MEM --only-error-correction" # --tmp-dir /tmp"

    J=$((${#F1[@]}-1))
    for I in `seq 0 $J`; do
        K=$((I+1))
        SPADESCMD=$SPADESCMD" --pe$K-1 ${F1[$I]} --pe$K-2 ${R1[$I]} "
    done

        echo -e $SPADESCMD

#	continue #for testing without running on cluster

        sbatch --partition=$NODE --ntasks=$CPUS --tasks-per-node=$CPUS --mem=$MEM"G" <<-EOF
#!/bin/bash
eval "$SPADESCMD"
EOF
    
done
