#/bin/bash

# this program merges prokka annotations and metabat depths
# make sure metabat $f1 dir has been renamed, runMetabat gives weird name

# input is one or more .fasta files
Bin="$@"

kaijudir=annotation/kaiju
prokkadir=annotation/prokka
metabatdir=annotation/metabat

for file in $Bin; do
        echo "start $file"

        f1=`basename $file .fasta`
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

# match orfs to nodes
awk -v OFS='\t' 'BEGIN{feature=0}{if($0 ~ /^>/){feature=$2};if($0 ~ /locus_tag/){print feature,$2}}' $prokkadir/prokka_$f1/$f1.tbl > $metabatdir/$f1.node_orf.tsv

#12000 missing contigs are cases where no CDS were called on the contig

# kaiju_nodes:	contigs:	178357
# kaiju_nodes:	orf:		433538

# now we want to merge prokka annotations with metabat depths

# add annotations from ffn
grep '^>' annotation/prokka/prokka_$f1/$f1.ffn | cut -f2 -d'>' | sed 's/ /'$'\t''/' > $metabatdir/$f1.ann.tsv
join -1 2 -2 1 -t$'\t' <(sort -k2,2 $metabatdir/$f1.node_orf.tsv) <(sort -k1,1 $metabatdir/$f1.ann.tsv) > $metabatdir/$f1.orf_node_ann.tsv

# merge metabat bin ids
grep -H '^>' $metabatdir/$f1/bin.*.fa | cut -f4 -d'/' | tr -s ':>' $'\t' > $metabatdir/$f1.binmap.tsv
join -1 2 -2 2 -t$'\t' -a2 -o 2.2,2.1,1.1,2.3 -e NA <(sort -k2,2 $metabatdir/$f1.binmap.tsv) <(sort -k2,2 $metabatdir/$f1.orf_node_ann.tsv) > $metabatdir/$f1.node_orf_bin_ann.tsv
# 433538

#get and sort metabat depths
join -1 1 -2 1 -t$'\t' <(sort -k1,1 -t$'\t' $metabatdir/$f1.node_orf_bin_ann.tsv) <(sort -k1,1 -t$'\t' $metabatdir/$f1.depth.txt) > $metabatdir/$f1.node_orf_bin_ann_depths.tsv
# 433535, missing 3 orfs
# 272917 hypothetical

done;


    
        
