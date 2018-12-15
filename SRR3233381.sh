#!/bin/bash
~/apps/bbmap/bbmap.sh in= in2= out=trimmed_SRR3233381_1.fastq.00.0_0.cor.fastq.gz.mapped.bam ref=./../spades-assembly/spades_SRR3233381_2018-10-20-13-56-31/scaffolds.fasta bamscript=bs.sh; sh bs.sh;
~/apps/bbmap/bbmap.sh in= in2= out=trimmed_SRR3233381_2.fastq.00.0_0.cor.fastq.gz.mapped.bam ref=./../spades-assembly/spades_SRR3233381_2018-10-20-13-56-31/scaffolds.fasta bamscript=bs.sh; sh bs.sh;
~/apps/bbmap/bbmap.sh in= in2= out=trimmed_SRR3233381__unpaired.00.0.mapped.bam ref=./../spades-assembly/spades_SRR3233381_2018-10-20-13-56-31/scaffolds.fasta bamscript=bs.sh; sh bs.sh;
~/apps/bbmap/bbmap.sh in= in2= out=trimmed_SRR3233383_1.fastq.00.0_0.cor.fastq.gz.mapped.bam ref=./../spades-assembly/spades_SRR3233381_2018-10-20-13-56-31/scaffolds.fasta bamscript=bs.sh; sh bs.sh;
~/apps/bbmap/bbmap.sh in= in2= out=trimmed_SRR3233383_2.fastq.00.0_0.cor.fastq.gz.mapped.bam ref=./../spades-assembly/spades_SRR3233381_2018-10-20-13-56-31/scaffolds.fasta bamscript=bs.sh; sh bs.sh;
~/apps/bbmap/bbmap.sh in= in2= out=trimmed_SRR3233383__unpaired.00.0.mapped.bam ref=./../spades-assembly/spades_SRR3233381_2018-10-20-13-56-31/scaffolds.fasta bamscript=bs.sh; sh bs.sh;
