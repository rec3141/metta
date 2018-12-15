#/bin/bash

# this program provides a CSV file that can be loaded into Bandage
# containing taxonomy labels and colors for each taxa
# DON'T RUN THIS IN PARALLEL on macpro, kaiju uses all the RAM

# input is one or more genome bins as .fasta files
Bin="$@"

kaijudir=annotation/kaiju
prokkadir=annotation/prokka
vizbindir=binning/vizbin
sketchdir=annotation/sketch
metabatdir=binning/metabat
spadesdir=assembly/spades

#pathway=sketch or pathway=prokka
pathway="prokka"

for file in $Bin; do
	echo "start $file"

	f1=`basename $file .scaffolds.fasta`
	echo $f1
	
	# #run minimap2/miniasm
#	/work/cryomics/apps/minimap2/minimap2 -X -t12 $file $file | gzip -1 > $f1.paf.gz
#	/work/cryomics/apps/miniasm/miniasm -c 1 -s 128 -h 100000 -f $file $f1.paf.gz > $f1.gfa

#	grep '^S'  $spadesdir/$f1/$f1.gfa | cut -f2,3 | sed 's/^\([0-9]\)/>\1/' | tr $'\t' $'\n' > $kaijudir/$f1.assembly.fasta

	cut -f1-3 $metabatdir/$f1.depth.txt | tr $'\t' ',' > $kaijudir/$f1.metabat.txt

#	ln -sf $file $f1.assembly.fasta
# to make prokka happy
	sed -E 's/(^>.{37}).*/\1/' $file > $kaijudir/$f1.assembly.fasta

    if [ "$pathway" == "sketch" ]; then
    #run sendsketch
    /work/cryomics/apps/bbmap/sendsketch.sh in=$file address=nt mode=sequence records=1 out=$sketchdir/$f1.sketch.tsv ow=t

	# match orfs to nodes
#	awk -v OFS=$'\t' 'BEGIN{feature=0}{if($0 ~ /^Query/){feature=$2};if($0 ~ /^[0-9]/){print feature,$12}}' $f1.sketch.tsv | sort | uniq -c | sed 's/^ *//' | sort -k2,2 -k1,1rn | uniq -f1 | sort -k2,2 | tr $'\t' ',' > $f1.nodes.csv
	awk -v OFS=$'\t' 'BEGIN{feature=0}{if($0 ~ /^Query/){feature=$2};if($0 ~ /^[0-9]/){print feature,$12 " " $13}}' $sketchdir/$f1.sketch.tsv | sort -k2,2 | tr $'\t' ',' > $sketchdir/$f1.nodes.csv

	#get orf colors
    cut -f2 -d',' $sketchdir/$f1.nodes.csv | uniq > $sketchdir/$f1.names.tsv
	color=`md5 <(echo -e $line) | cut -f4 -d' '| cut -b1-6`
	while read line; do echo "#$color";  done < $sketchdir/$f1.names.tsv > $sketchdir/$f1.colors.tsv
	paste -d',' $sketchdir/$f1.colors.tsv $sketchdir/$f1.names.tsv > $sketchdir/$f1.taxcolors.csv
	
	#merge orfs and colors
	echo -e "Node,Taxonomy,Color" > $sketchdir/$f1.bandage.csv
    join -1 2 -2 2 -t',' $sketchdir/$f1.nodes.csv $sketchdir/$f1.taxcolors.csv | awk -v OFS=',' -v FS=',' '{print $2,$1,$3}' >> $sketchdir/$f1.bandage.csv
  
    #merge with metabat
    sort -k1,1 -t',' $kaijudir/$f1.metabat.txt > $sketchdir/tmp1
    sort -k1,1 -t',' $sketchdir/$f1.bandage.csv > $sketchdir/tmp2
    join -1 1 -2 1 -t',' $sketchdir/tmp1 $sketchdir/tmp2 > $sketchdir/$f1.rcolor.csv

    fi;

    if [ "$pathway" == "prokka" ]; then

	#run prokka
	if [ ! -d "$prokkadir/prokka_$f1" ]; then 
		echo "running prokka"
		/work/cryomics/apps/prokka-biolinux/prokka/bin/prokka --metagenome --outdir $prokkadir/prokka_$f1 --prefix $f1 --locustag $f1 --force --fast --rawproduct $kaijudir/$f1.assembly.fasta
	fi;

	# match orfs to nodes
	awk -v OFS='\t' 'BEGIN{feature=0}{if($0 ~ /^>/){feature=$2};if($0 ~ /locus_tag/){print feature,$2}}' $prokkadir/prokka_$f1/$f1.tbl | sort -k2,2 > $kaijudir/$f1.kaiju_nodes.tsv
        
	#run kaiju
	if [ ! -e "$kaijudir/$f1.kaiju.tsv" ]; then
		echo "running kaiju"
		/work/cryomics/apps/kaiju-1.6.3/bin/kaiju -t /scratch/cryomics/reference_dbs/kaiju/2017/nodes.dmp -f /scratch/cryomics/reference_dbs/kaiju/2017/kaiju_db_nr_euk.fmi -p -i $prokkadir/prokka_$f1/$f1.faa -a greedy -e 5 -o $kaijudir/$f1.kaiju.tsv -z 12
	fi;
	
	#get orf taxonomy
   	/work/cryomics/apps/kaiju-1.6.3/bin/addTaxonNames -t /scratch/cryomics/reference_dbs/kaiju/2017/nodes.dmp -n /scratch/cryomics/reference_dbs/kaiju/2017/names.dmp -r "no rank",superkingdom,phylum,class,order,family,genus,species -i $kaijudir/$f1.kaiju.tsv -o $kaijudir/$f1.kaiju_names_all.tsv
# options: class cohort family forma genus infraclass infraorder kingdom no rank order parvorder phylum species species group species subgroup subclass subfamily subgenus subkingdom suborder subphylum subspecies subtribe superclass superfamily superkingdom superorder superphylum tribe varietas

    #sort metabat
    sort -k1,1 -t',' $kaijudir/$f1.metabat.txt > $kaijudir/$f1.metabat.sort.txt

    #taxonomy list
    for tax in superkingdom phylum class order family genus species; do

        echo $tax
    	/work/cryomics/apps/kaiju-1.6.3/bin/addTaxonNames -t /scratch/cryomics/reference_dbs/kaiju/2017/nodes.dmp -n /scratch/cryomics/reference_dbs/kaiju/2017/names.dmp -r $tax -i $kaijudir/$f1.kaiju.tsv -o $kaijudir/$f1.kaiju_names_$tax.tsv

        sort -k3 $kaijudir/$f1.kaiju_names_$tax.tsv | cut -f3,4 | sort -k1,1 -u | tr $'\t' '|' > $kaijudir/$f1.kaiju_names_sort_$tax.tsv

        #merge orfs and nodes, keeping only most common taxonomic hit
        join -1 2 -2 2 -t $'\t' $kaijudir/$f1.kaiju_nodes.tsv <(sort -k2,2 $kaijudir/$f1.kaiju_names_$tax.tsv) | cut -f2,4 | sort | uniq -c | sed 's/^ *//' | sort -k2,2 -k1,1rn | awk '{print $3 "\t" $2}' |  grep -v '^0' | uniq -f1 | sort -k1,1 | tr $'\t' '|' > $kaijudir/$f1.nodes_$tax.tsv
        join -1 1 -2 1 -t '|' <(sort -k1,1 $kaijudir/$f1.nodes_$tax.tsv) <(sort -k1,1 $kaijudir/$f1.kaiju_names_sort_$tax.tsv) | cut -f2- -d'|'| tr '|' ',' | sort -k2 -t',' >  $kaijudir/$f1.nodes_$tax.csv


        #get orf colors
        cut -f2 -d',' $kaijudir/$f1.nodes_$tax.csv | uniq > $kaijudir/$f1.names_$tax.tsv
		color=`md5 <(echo -e $line) | cut -f4 -d' '| cut -b1,3,5,7,9,11`
        while read line; do echo "#$color";  done < $kaijudir/$f1.names_$tax.tsv > $kaijudir/$f1.colors_$tax.tsv
        paste -d',' $kaijudir/$f1.colors_$tax.tsv $kaijudir/$f1.names_$tax.tsv > $kaijudir/$f1.taxcolors_$tax.csv
    
        #merge orfs and colors
        echo -e "Node,Taxonomy,Color" > $kaijudir/$f1.bandage_$tax.csv
        join -1 2 -2 2 -t',' <(sort -k2,2 -t',' $kaijudir/$f1.nodes_$tax.csv) <(sort -k2,2 -t',' $kaijudir/$f1.taxcolors_$tax.csv) | awk -v OFS=',' -v FS=',' '{print $2,$1,$3}' >> $kaijudir/$f1.bandage_$tax.csv

        #merge with metabat
        sort -k1,1 -t',' $kaijudir/$f1.bandage_$tax.csv > $kaijudir/$f1.bandage_$tax.sort.csv
        join -1 1 -2 1 -t',' $kaijudir/$f1.metabat.sort.txt $kaijudir/$f1.bandage_$tax.sort.csv > $kaijudir/$f1.rcolor_$tax.csv

    done;

    fi;
   
   	echo "done $file"
done;




