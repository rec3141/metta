#!/bin/bash
# this program gets the reads that map to particular contigs
# input $1 is assembly with contigs to get reads for
# input $2 is directory containing bam files

assembly=$1
f1=`basename $assembly .fasta`
bamdir=$2
mappeddir=mapped/$f1
mkdir -p $mappeddir

for bamfile in $bamdir/*.bam; do
b1=`basename $bamfile .mapped_sorted.bam`

  #samtools view -b trimmed_Plate1-11A-18_S81.mapped_sorted.bam NODE_2_length_478704_cov_10.545229 | samtools fastq -1 r1.fastq.gz -2 r2.fastq.gz -s r0.fastq.gz -
  CMD="samtools view -b $bamfile \`cat <(grep '>' $assembly | cut -f2 -d'>' | tr $'\n' ' ')\` | samtools fastq -1 $mappeddir/$f1\"_\"$b1\"_R1.fastq.gz\" -2 $mappeddir/$f1\"_\"$b1\"_R2.fastq.gz\" -s $mappeddir/$f1\"_\"$b1\"_R0.fastq.gz\" -"
  echo $CMD

#bbsplit.sh in1=reads1.fq in2=reads2.fq ref=ecoli.fa,salmonella.fa basename=out_%.fq outu1=clean1.fq outu2=clean2.fq

done


#~/apps/SPAdes-3.13.0-Linux/bin/spades.py -o assembly/spades/spades_cyanos -t 7 -m 400 --tmp-dir /tmp --only-assembler --careful -1 mapped/cyanos/cyanos_R1.fastq.gz -2 mapped/cyanos/cyanos_R2.fastq.gz -s mapped/cyanos/cyanos_R0.fastq.gz

