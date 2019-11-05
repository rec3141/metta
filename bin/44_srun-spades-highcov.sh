#!/bin/bash
#input $1 is directory with corrected reads
NODE=bio
THREADS=7
MEM=500

corrected=$1

#trimmed_Plate1-1F-15_S6_L001_R1_001.fastq.00.0_0.cor.fastq.gz
r1=`basename $1`

    OUTDIR=all_highcov_$r1"_"`date +'%F-%H-%M-%S'`

    SPADESCMD="~/apps/SPAdes-3.13.0-Linux/bin/spades.py -o assembly/spades/spades_$OUTDIR -t $THREADS -m $MEM --tmp-dir /tmp --only-assembler --12 $corrected/all_pairs_highcov.fq.gz --merged $corrected/all_merged_highcov.fq.gz -s $corrected/all_singles_highcov.fq.gz"

    echo -e $SPADESCMD

        sbatch --partition=$NODE --ntasks=$THREADS --tasks-per-node=$THREADS --mem=$MEM"G" <<-EOF
#!/bin/bash
eval "$SPADESCMD"
EOF
