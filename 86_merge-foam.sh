#/bin/bash
# input 1 is nodes file; e.g. annotation/kaiju/scaffolds.1000bpp.kaiju_nodes.tsv
# input 2 is depth file; e.g. binning/metabat/scaffolds.fasta.depth.txt

kaijudir="annotation/kaiju"
foamdir="annotation/foam"

nodefile=$1
depthfile=$2

f1=`basename $nodefile .kaiju_nodes.tsv`

#foam
join -1 2 -2 1 -t $'\t' <(sort -k2,2 -t $'\t' $nodefile) <(sort -k1,1 -t $'\t' $foamdir/$f1.foam.tsv) > $foamdir/$f1.foam.ont.tsv
join -1 2 -2 1 -t $'\t' <(sort -k2,2 -t $'\t' $foamdir/$f1.foam.ont.tsv) <(sort -k1,1 -t $'\t' $depthfile) > $kaijudir/$f1.foam.depths.tsv
echo "contig orf score foam_hmm ko foam_1 foam_2 foam_3 foam_4 "`head -n1 $depthfile | cut -f2-` | sed 's/trimmed_//g;s/.mapped_sorted.bam//g' > $foamdir/foam.headers.txt

#kegg
join -1 2 -2 1 -t $'\t' <(sort -k2,2 -t $'\t' $nodefile) <(sort -k1,1 -t $'\t' $foamdir/$f1.kegg.tsv) > $foamdir/$f1.kegg.ont.tsv
join -1 2 -2 1 -t $'\t' <(sort -k2,2 -t $'\t' $foamdir/$f1.kegg.ont.tsv) <(sort -k1,1 -t $'\t' $depthfile) > $kaijudir/$f1.kegg.depths.tsv
echo "contig orf score foam_hmm ko ko_1 ko_2 ko_3 gene description ec "`head -n1 $depthfile | cut -f2-` | sed 's/trimmed_//g;s/.mapped_sorted.bam//g' > $foamdir/kegg.headers.txt


#join -1 2 -2 1 -t $'\t' <(sort -k2,2 -t $'\t' $nodefile) <(sort -k1,1 -t $'\t' $ontfile) > $foamdir/$f2.ont.tsv
#join -1 2 -2 1 -t $'\t' <(sort -k2,2 -t $'\t' $foamdir/$f2.ont.tsv) <(sort -k1,1 -t $'\t' $depthfile) > $kaijudir/$f2.depths.tsv
