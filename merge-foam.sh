#/bin/bash
# input 1 is nodes file; e.g. annotation/kaiju/scaffolds.1000bpp.kaiju_nodes.tsv
# input 2 is ontology file; e.g.  annotation/foam/scaffolds.1000bpp.kegg.tsv
# input 3 is depth file; e.g. binning/metabat/scaffolds.fasta.depth.txt

kaijudir="annotation/kaiju"
foamdir="annotation/foam"

nodefile=$1
ontfile=$2
depthfile=$3

f1=`basename $nodefile .kaiju_nodes.tsv`
f2=`basename $ontfile .tsv`

join -1 2 -2 1 -t $'\t' <(sort -k2,2 -t $'\t' $nodefile) <(sort -k1,1 -t $'\t' $ontfile) > $foamdir/$f2.ont.tsv

join -1 2 -2 1 -t $'\t' <(sort -k2,2 -t $'\t' $foamdir/$f2.ont.tsv) <(sort -k1,1 -t $'\t' $depthfile) > $kaijudir/$f2.depths.tsv
