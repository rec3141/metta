#!/bin/bash

# this program parses the output from hmmer 
# input 1 is .faa sequence file

SEQFILE=$1

for INDIR in `ls -d */ | cut -f1 -d'/'`; do

	rm $INDIR/$INDIR.$sf.hits.txt

	sf=`basename $SEQFILE .faa`

	for FILE in $INDIR/tsv/*.$sf.tsv; do 
		echo $FILE; 
		cat $FILE | grep -v '^#' | sed -E 's/ +/\t/g;' | cut -f1,4,6 | sort -k1,1 -k3,3rn | awk -v OFS=$'\t' '{print $3,$2,$1}' | uniq -f2 >> $INDIR/$INDIR.$sf.hits.txt
	done

done

cat */*.$sf.hits.txt | sort -k3,3 -k1,1rn | uniq -f2 | sort -k2,2 > $sf.hits.txt

#merge with KEGG ontology
join -1 1 -2 2 -a 2 -e NA -t $'\t' -o 2.3,2.1,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8 KEGG-hmm-ont.txt $sf.hits.txt | sort -u | sort -k1,1 > $sf.kegg.tsv

#merge with FOAM ontology
join -1 1 -2 2 -a 2 -e NA -t $'\t' -o 2.3,2.1,2.2,1.2,1.3,1.4,1.5,1.6  FOAM-hmm-ont.txt $sf.hits.txt | sort -u | sort -k1,1 > $sf.foam.tsv


# this cuts out the sequence, profile, and score
# then sort them and outputs only the best profile hit for each sequence
# 386.3 HMMsoil25897 assembly_graph_ns_k127_00032
# 180.7 HMMsoil24834 assembly_graph_ns_k127_00052
# 56.7 HMMsoil26923 assembly_graph_ns_k127_00054
# 18.2 HMMsoil5120 assembly_graph_ns_k127_00072

# merge with prokka output to get which contig each hmm hit is on, then get taxonomy
# contig	orf	taxonomy	bin?	hmmhit	score	FOAM1	FOAM2	FOAM3	FOAM4	FOAM5	sample1-cov	sample2-cov	sample3-cov	...
