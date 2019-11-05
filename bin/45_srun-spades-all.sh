#!/bin/bash
#input $1 is directory of corrected reads
NODE=bio
THREADS=7
MEM=1500

corrected=$1

#trimmed_Plate1-1F-15_S6_L001_R1_001.fastq.00.0_0.cor.fastq.gz
cat $corrected/*_R1_*.gz $corrected/*_1.*.gz > $corrected/R1.cor.fastq.gz
cat $corrected/*_R2_*.gz $corrected/*_2.*.gz > $corrected/R2.cor.fastq.gz
cat $corrected/*_R_.*.gz $corrected/*__.*.gz > $corrected/R0.cor.fastq.gz

    F1=$corrected/R1.cor.fastq.gz
    R1=$corrected/R2.cor.fastq.gz
    S1=$corrected/R0.cor.fastq.gz

    OUTDIR=all"_"`date +'%F-%H-%M-%S'`

    SPADESCMD="~/apps/SPAdes-3.13.0-Linux/bin/spades.py -o assembly/spades/spades_$OUTDIR -t $THREADS -m $MEM --tmp-dir /tmp --only-assembler -1 $F1 -2 $R1 -s $S1"

    echo -e $SPADESCMD

        sbatch --partition=$NODE --ntasks=$THREADS --tasks-per-node=$THREADS --mem=$MEM"G" <<-EOF
#!/bin/bash
eval "$SPADESCMD"
EOF
