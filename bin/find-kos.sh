#!/bin/bash
#for dir in `cut -f1 FOAM-onto_rel1.tsv | cut -f1 -d' ' | uniq`; do ./find-kos.sh $dir; done

Search=$1

grep -i "$Search" FOAM-onto_rel1.tsv | sort -u > $Search.ont
cut -f5 $Search.ont > $Search.ko
wc -l $Search.ko
grep -A1 -f $Search.ko FOAM-hmm_rel1a.hmm | grep NAME | tr -s ' ' | cut -f2 -d' ' | sort -u > $Search.acc
wc -l $Search.acc
~/apps/hmmer/binaries/hmmfetch -f FOAM-hmm_rel1a.hmm $Search.acc > tmp.$Search
hmms=`grep -c ACC tmp.$Search`
pieces=`expr $hmms / 1000 + 1`
echo $pieces

 for i in `seq 1 $pieces`; do
	echo $i
	awk -v i=$i 'BEGIN{RS="//";ORS="//";start=1000*i-1000;end=1000*i} {if(NR > start && NR <= end) print $0}' tmp.$Search | uniq > $Search.$i.hmm
 done;

#should remove trailing // on final file to avoid crashes (uniq)
