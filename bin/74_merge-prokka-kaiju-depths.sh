#/bin/bash
# this program merges prokka annotations, kaiju taxonomy, and metabat depths

# e.g. bin/74_merge-prokka-kaiju-depths.sh assembly/spades/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/scaffolds_gt2000.fasta annotation/metabat/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/ASGARD_20181114/scaffolds_gt2000.fasta.depth.txt
# input $1 is the metabat depth matrix
# input $2 is directory containing metabat bins 

kaijudir=annotation/kaiju
prokkadir=annotation/prokka
kjdb=~/work/sfosdna/reference_dbs/kaiju
kaijubin=~/apps/kaiju/bin

depthfile=$1
bindir=$2

f1=`basename $depthfile .fasta.depth.txt`
echo $f1

echo "a"
#get orf taxonomy
$kaijubin/addTaxonNames -t $kjdb/nodes.dmp -n $kjdb/names.dmp -r "no rank",superkingdom,phylum,class,order,family,genus,species -i $kaijudir/$f1.kaiju.tsv -o $kaijudir/$f1.kaiju_names_all.tsv
# options: class cohort family forma genus infraclass infraorder kingdom no rank order parvorder phylum species species group species subgroup subclass subfamily subgenus subkingdom suborder subphylum subspecies subtribe superclass superfamily superkingdom superorder superphylum tribe varietas

# match orfs to nodes
awk -v OFS='\t' 'BEGIN{feature=0}{if($0 ~ /^>/){feature=$2};if($0 ~ /locus_tag/){print feature,$2}}' $prokkadir/prokka_$f1/$f1.tbl | sort -k2,2 > $kaijudir/$f1.kaiju_nodes.tsv

echo "b"

#sort depth matrix
cut -f1-3 $depthfile | tr $'\t' ',' > $kaijudir/$f1.meandepth.txt
sort -k1,1 -t',' $kaijudir/$f1.meandepth.txt > $kaijudir/$f1.meandepth.sort.txt

echo "c"

#taxonomy list
for tax in superkingdom phylum class order family genus species; do

	echo $tax
	$kaijubin/addTaxonNames -t $kjdb/nodes.dmp -n $kjdb/names.dmp -r $tax -i $kaijudir/$f1.kaiju.tsv -o $kaijudir/$f1.kaiju_names_$tax.tsv

	# 3 is ncbid taxonomy; 4 is score
#	sort -k3 $kaijudir/$f1.kaiju_names_$tax.tsv | cut -f3,4 | sort -k1,1 -u | tr $'\t' '|' > $kaijudir/$f1.kaiju_names_sort_$tax.tsv
	sort -k3 $kaijudir/$f1.kaiju_names_$tax.tsv | cut -f3,8 | sort -k1,1 -u | tr $'\t' '|' > $kaijudir/$f1.kaiju_names_sort_$tax.tsv

echo "d"

	#merge orfs and nodes, keeping only most common taxonomic hit
	join -1 2 -2 2 -t $'\t' <(sort -k2,2 -t $'\t' $kaijudir/$f1.kaiju_nodes.tsv) <(sort -k2,2 -t $'\t' $kaijudir/$f1.kaiju_names_$tax.tsv) | cut -f2,4 | sort | uniq -c | sed 's/^ *//' | sort -k2,2 -k1,1rn | awk '{print $3 "\t" $2}' |  grep -v '^0' | uniq -f1 | sort -k1,1 | tr $'\t' '|' > $kaijudir/$f1.nodes_$tax.tsv
	join -1 1 -2 1 -t '|' <(sort -k1,1 -t'|' $kaijudir/$f1.nodes_$tax.tsv) <(sort -k1,1 -t '|' $kaijudir/$f1.kaiju_names_sort_$tax.tsv) | cut -f2- -d'|' | tr '|' ',' | sort -k2 -t',' >  $kaijudir/$f1.nodes_$tax.csv

echo "e"

	#get orf colors
	cut -f2 -d',' $kaijudir/$f1.nodes_$tax.csv | uniq > $kaijudir/$f1.names_$tax.tsv
	while read line; do echo "#"`md5sum <(echo -e $line) | cut -f1 -d' '| cut -b1,3,5,7,9,11`;  done < $kaijudir/$f1.names_$tax.tsv > $kaijudir/$f1.colors_$tax.tsv
	paste -d',' $kaijudir/$f1.colors_$tax.tsv $kaijudir/$f1.names_$tax.tsv > $kaijudir/$f1.taxcolors_$tax.csv

echo "f"

	#merge orfs and colors
	echo -e "Node,Taxonomy,Color" > $kaijudir/$f1.bandage_$tax.csv
	join -1 2 -2 2 -t',' <(sort -k2,2 -t',' $kaijudir/$f1.nodes_$tax.csv) <(sort -k2,2 -t',' $kaijudir/$f1.taxcolors_$tax.csv) | awk -v OFS=',' -v FS=',' '{print $2,$1,$3}' >> $kaijudir/$f1.bandage_$tax.csv

echo "g"

	#merge with metabat
	sort -k1,1 -t',' $kaijudir/$f1.bandage_$tax.csv > $kaijudir/$f1.bandage_$tax.sort.csv
	join -1 1 -2 1 -t',' <(sort -k1,1 -t',' $kaijudir/$f1.meandepth.sort.txt) <(sort -k1,1 -t',' $kaijudir/$f1.bandage_$tax.sort.csv) > $kaijudir/$f1.rcolor_$tax.csv

done;



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
	grep -H '^>' $binfile | rev | cut -f1 -d'/' | rev | tr -s ':>' $'\t' >> $kaijudir/$f1.binmap.tsv
done < <( find $bindir -name "$f1*.fa" )

join -1 2 -2 2 -t$'\t' -a2 -o 2.2,2.1,1.1,2.3 -e NA <(sort -k2,2 -t$'\t' $kaijudir/$f1.binmap.tsv) <(sort -k2,2 -t$'\t' $kaijudir/$f1.orf_node_ann.tsv) > $kaijudir/$f1.node_orf_bin_ann.tsv
# 433538
# some contigs showing up NA for bins but all contigs should have been assigned a bin... 

echo "c"
#get and sort metabat depths
join -1 1 -2 1 -t$'\t' <(sort -k1,1 -t$'\t' $kaijudir/$f1.node_orf_bin_ann.tsv) <(sort -k1,1 -t$'\t' $depthfile) > $kaijudir/$f1.node_orf_bin_ann_depths.tsv
# 433535, missing 3 orfs
# 272917 hypothetical

echo "d"
# add taxonomy
join -1 1 -2 2 -t$'\t' <(cat $kaijudir/$f1.kaiju_names_all.tsv | cut -f2,8 | awk 'BEGIN{FS="\t";OFS="\t"} {if(NF==1) {print $0,"unclassified;"} else print $0}' | sort -k1,1 -t$'\t') <(sort -k2,2 -t$'\t' $kaijudir/$f1.node_orf_bin_ann_depths.tsv) > $kaijudir/$f1.node_orf_bin_ann_depths_tax.tsv


echo "done $f1"




# contigs missing bins
# NODE_100786_length_2233_cov_1.786501
# NODE_101239_length_2226_cov_2.258406
# NODE_101266_length_2225_cov_5.711982
# NODE_101472_length_2222_cov_2.140286
# NODE_105672_length_2156_cov_2.853879
# NODE_112486_length_2059_cov_2.626248
# NODE_113200_length_2049_cov_2.687563
# NODE_113292_length_2048_cov_2.002007
# NODE_113295_length_2048_cov_1.701455
# NODE_113437_length_2046_cov_2.326971
# NODE_11524_length_11376_cov_38.008656
# NODE_116051_length_2012_cov_3.125703
# NODE_12240_length_10834_cov_41.551999
# NODE_14764_length_9433_cov_92.395607
# NODE_1846_length_40140_cov_83.226868
# NODE_22182_length_6939_cov_12.927804
# NODE_23818_length_6579_cov_85.716891
# NODE_24831_length_6379_cov_3.196237
# NODE_2504_length_33538_cov_42.809874
# NODE_4024_length_24600_cov_77.539988
# NODE_4530_length_22531_cov_84.628359
# NODE_47465_length_3924_cov_5.449212
# NODE_47825_length_3903_cov_38.159823
# NODE_50854_length_3722_cov_16.275430
# NODE_5680_length_19170_cov_95.130055
# NODE_578_length_75127_cov_40.816656
# NODE_601_length_73822_cov_14.749251
# NODE_61515_length_3226_cov_16.898770
# NODE_64984_length_3099_cov_2.315046
# NODE_66213_length_3054_cov_43.621541
# NODE_74081_length_2808_cov_2.752997
# NODE_7546_length_15577_cov_24.255444
# NODE_75572_length_2768_cov_2.414670
# NODE_79539_length_2662_cov_2.422325
# NODE_79848_length_2654_cov_2.167757
# NODE_81123_length_2624_cov_2.595952
# NODE_81340_length_2619_cov_1.813183
# NODE_86073_length_2509_cov_3.033007
# NODE_8765_length_13936_cov_79.585044
# NODE_92281_length_2385_cov_2.195708
# NODE_94608_length_2342_cov_2.192392
# NODE_95737_length_2321_cov_2.452339
# NODE_99056_length_2262_cov_1.846851




