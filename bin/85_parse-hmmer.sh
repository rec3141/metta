#!/bin/bash

# this program parses the output from hmmer 
# input $1 is .faa sequence file

seqfile=$1
foamdir=annotation/foam

for hmmd in `ls -d $foamdir/*/ `; do

	h1=`basename $hmmd`
	s1=`basename $seqfile .faa`

	rm $hmmd/$h1.$s1.hits.txt

	for tsvfile in $hmmd/tsv/*.$s1.tsv; do 
		echo $tsvfile;	
		cat $tsvfile | grep -v '^#' | sed -E 's/ +/\t/g;' | cut -f1,4,6 | sort -k1,1 -k3,3rn | awk -v OFS=$'\t' '{print $3,$2,$1}' | uniq -f2 >> $hmmd/$h1.$s1.hits.txt
	done
done

cat $foamdir/*/*.$s1.hits.txt | sort -k3,3 -k1,1rn | uniq -f2 | sort -k2,2 > $foamdir/$s1.hits.txt

#merge with KEGG ontology
join -1 1 -2 2 -a 2 -e NA -t $'\t' -o 2.3,2.1,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8 $foamdir/KEGG-hmm-ont.txt $foamdir/$s1.hits.txt | sort -u | sort -k1,1 > $foamdir/$s1.kegg.tsv

#merge with FOAM ontology
join -1 1 -2 2 -a 2 -e NA -t $'\t' -o 2.3,2.1,2.2,1.2,1.3,1.4,1.5,1.6  $foamdir/FOAM-hmm-ont.txt $foamdir/$s1.hits.txt | sort -u | sort -k1,1 > $foamdir/$s1.foam.tsv


# this cuts out the sequence, profile, and score
# then sort them and outputs only the best profile hit for each sequence
# 386.3 HMMsoil25897 assembly_graph_ns_k127_00032
# 180.7 HMMsoil24834 assembly_graph_ns_k127_00052
# 56.7 HMMsoil26923 assembly_graph_ns_k127_00054
# 18.2 HMMsoil5120 assembly_graph_ns_k127_00072

# merge with prokka output to get which contig each hmm hit is on, then get taxonomy
# contig	orf	taxonomy	bin?	hmmhit	score	FOAM1	FOAM2	FOAM3	FOAM4	FOAM5	sample1-cov	sample2-cov	sample3-cov	...
