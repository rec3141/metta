#!/bin/bash
#input $1 is assembly fasta file
#input $2 is directory containing sorted bam files

assembly=$1
bamdir=$2

sbatch --partition=bio --ntasks=7 --tasks-per-node=7 --mem=1400G <<-EOF
#!/bin/bash
#~/apps/metabat/runMetaBat.sh -m 2000 --maxP 99 --minS 2 -s 2000 --maxEdges 1000 --minCV 0 --saveCls -t7 -v -d --unbinned $assembly $bamdir/*sorted.bam
#~/apps/metabat/metabat2  -m 2000 --maxP 99 --minS 2 -s 2000 --maxEdges 1000 --minCV 0 --saveCls -t7 -v -d --unbinned --inFile assembly/spades/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/scaffolds_gt2000.fasta --outFile scaffolds_gt2000.fasta.metabat-bins --unbinned/bin --abdFile scaffolds_gt2000.fasta.depth.txt
~/apps/metabat/metabat2  -m 2000 --maxP 95 --minS 60 -s 2000 --maxEdges 200 --minCV 1 --saveCls -t7 -v -d --unbinned --inFile assembly/spades/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/scaffolds_gt2000.fasta --outFile annotation/metabat/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/ASGARD_20181114/`date +'%F-%H-%M-%S'`/scaffolds_gt2000.fasta.metabat-bins --unbinned/bin --abdFile annotation/metabat/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/ASGARD_20181114/scaffolds_gt2000.fasta.depth.txt
#mkdir -p annotation/metabat/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/ASGARD_20181114
#mv scaff* annotation/metabat/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/ASGARD_20181114/
EOF


