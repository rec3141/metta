#/bin/bash

# this program merges prokka annotations and metabat depths
#e.g. bin/74_merge-prokka-kaiju-depths.sh assembly/spades/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/scaffolds_gt2000.fasta annotation/metabat/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/ASGARD_20181114/2018-12-07-16-32-16

# input $1 is the metabat depth matrix
# input $2 is directory containing metabat bins 

kaijudir=annotation/kaiju
prokkadir=annotation/prokka

depthfile=$1
metabatdir=$2

f1=`basename $depthfile .fasta.depth.txt`
echo $f1

# contigs:	spades:	190644
# 
# contigs:	fna:	190644
# CDS:		ffn:	433538
# proteins:	faa:	428693
# 
# grep -c 'locus_tag' annotation/prokka/prokka_scaffolds.1000bpp/scaffolds.1000bpp.tbl
# grep '^[0-9]' annotation/prokka/prokka_scaffolds.1000bpp/scaffolds.1000bpp.tbl | cut -f3 | grep -c 'repeat_region'
# 
# contigs:	tbl:	190644
# CDS:		tbl:	428693
# rRNA:		tbl:	385
# tRNA:		tbl:	4366
# tmRNA:		tbl:	94
# 
# not counted
# repeat_region		36 (CRISPR)
# 
# cut -f1 annotation/kaiju/scaffolds.1000bpp.kaiju_nodes.tsv  | sort -u | wc -l
# cut -f2 annotation/kaiju/scaffolds.1000bpp.kaiju_nodes.tsv  | sort -u | wc -l
#12000 missing contigs are cases where no CDS were called on the contig

# kaiju_nodes:	contigs:	178357
# kaiju_nodes:	orf:		433538

# now we want to merge prokka annotations with metabat depths

echo "a"
# add annotations from ffn
grep '^>' $prokkadir/prokka_$f1/$f1.ffn | cut -f2 -d'>' | sed 's/ /'$'\t''/' > $kaijudir/$f1.ann.tsv
join -1 2 -2 1 -t$'\t' <(sort -k2,2 $kaijudir/$f1.kaiju_nodes.tsv) <(sort -k1,1 $kaijudir/$f1.ann.tsv) > $kaijudir/$f1.orf_node_ann.tsv

echo "b"
# merge metabat bin ids
# scaffolds_gt2000.fasta.metabat-bins.1002.fa
rm $kaijudir/$f1.binmap.tsv
while read binfile; do
	grep -H '^>' $binfile | xargs -I{} basename {} | tr -s ':>' $'\t' >> $kaijudir/$f1.binmap.tsv
done < <( find $metabatdir -name "$f1*.fa" )

join -1 2 -2 2 -t$'\t' -a2 -o 2.2,2.1,1.1,2.3 -e NA <(sort -k2,2 -t$'\t' $kaijudir/$f1.binmap.tsv) <(sort -k2,2 -t$'\t' $kaijudir/$f1.orf_node_ann.tsv) > $kaijudir/$f1.node_orf_bin_ann.tsv
# 433538

echo "c"
#get and sort metabat depths
join -1 1 -2 1 -t$'\t' <(sort -k1,1 -t$'\t' $kaijudir/$f1.node_orf_bin_ann.tsv) <(sort -k1,1 -t$'\t' $depthfile) > $kaijudir/$f1.node_orf_bin_ann_depths.tsv
# 433535, missing 3 orfs
# 272917 hypothetical


    
        
