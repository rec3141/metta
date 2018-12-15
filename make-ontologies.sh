#download KEGG Ontology database
#wget "https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir=" -O KEGG.txt
#wget "https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir=" -O KEGG.json

#clean up
grep name KEGG.json | sed -E 's/'$'\t''/#/g;s/"//g;s/name://g;s/  /'$'\t''/g;s/; /'$'\t''/g;s/ \[EC:(.*)\]/'$'\t''\1/g;s#\\/#/#g;s/,$//;s/'\''//g;' > KEGG.tsv

#make full KEGG ontology
awk -F "[#'$'\t'']" -v OFS=$'\t' 'BEGIN{A=0;B=0;C=0;D=0;};
    {
    if($0 ~ "^##0") {A=$3};
    if($0 ~ "^###0") {B=$4};
    if($0 ~ "^####0") {C=$5};
    if($0 ~ "^#####K") {print $6,A,B,C,$7,$8,$9}
    }' KEGG.tsv | sort -u > KEGG-ont.tsv

#extract KOs from FOAM
grep -B1 ACC FOAM-hmm_rel1a.hmm | paste - - - | tr -s ' ' | tr ' ' $'\t' | cut -f2,4 | awk '{print $2,$1}' | cut -f1 -d'_' | sed 's/KO://g' | awk 'BEGIN{ FS="[ ,]"; OFS="\t" }; { for (i=2;i<=NF;i+=1) {printf "%s%s%s%s", $1, OFS, $i, ORS}}' | sort -k2,2 > FOAM-hmm-ko.txt

#merge FOAM hmms and ontology
join -1 5 -2 2 -t $'\t' -o 2.1,2.2,1.1,1.2,1.3,1.4 -e NA  <(sort -k5,5 -t $'\t' FOAM-onto_rel1.tsv) FOAM-hmm-ko.txt | sort -k1,1 | uniq  > FOAM-hmm-ont.txt

#merge FOAM and KEGG
join -1 1 -2 2 -t $'\t' -a2 -e NA -o 2.1,2.2,1.2,1.3,1.4,1.5,1.6,1.7 KEGG-ont.tsv FOAM-hmm-ko.txt | sort -k1,1 | uniq > KEGG-hmm-ont.txt

rm FOAM-hmm-ko.txt
rm KEGG-ont.tsv
