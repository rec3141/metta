#!/bin/bash
#input $1 is directory of corrected reads
#input $2 is mindepth level 0 1 2 3

NODE=bio
THREADS=7
MEM=1400

corrected=$1
mindepth=$2

    pairs=$corrected/norm_$mindepth"_"pairs.fastq.gz
    merged=$corrected/norm_$mindepth"_"merged.fastq.gz
    singles=$corrected/norm_$mindepth"_"singles.fastq.gz

    c1=`basename $corrected`

    OUTDIR="all_norm_"$c1"_"`date +'%F-%H-%M-%S'`

#    SPADESCMD="~/apps/SPAdes-3.13.0-Linux/bin/spades.py -o assembly/spades/spades_$OUTDIR -t $THREADS -m $MEM --tmp-dir /tmp --only-assembler --careful --12 $pairs -s $singles --merged $merged"
    SPADESCMD="~/apps/SPAdes-3.13.0-Linux/bin/spades.py --meta -o assembly/spades/spades_$OUTDIR -t $THREADS -m $MEM --tmp-dir /tmp --only-assembler --12 $pairs --merged $merged"

    echo -e $SPADESCMD

        sbatch --partition=$NODE --ntasks=$THREADS --tasks-per-node=$THREADS --mem=$MEM"G" <<-EOF
#!/bin/bash
eval "$SPADESCMD"
EOF
