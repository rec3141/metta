#!/bin/bash
NODE=bio
THREADS=7
MEM=314

    F1=corrected/R1.cor.fastq.gz
    R1=corrected/R2.cor.fastq.gz
    S1=corrected/R0.cor.fastq.gz

    OUTDIR=all"_"`date +'%F-%H-%M-%S'`

#    SPADESCMD="~/apps/SPAdes-3.12.0-Linux/bin/spades.py -o spades-assembly/spades_$OUTDIR -t $THREADS -m $MEM --tmp-dir /tmp --only-assembler -1 $F1 -2 $R1 -s $S1"
    PLASSCMD="plass assemble $F1 $R1 plass-assembly/$OUTDIR/plass-assembly.faa tmp"
    echo -e $SPADESCMD

        sbatch --partition=$NODE --ntasks=$THREADS --tasks-per-node=$THREADS --mem=$MEM"G" <<-EOF
#!/bin/bash
eval "$PLASSCMD"
EOF
