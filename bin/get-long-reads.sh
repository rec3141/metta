perl -0076 -ne 'chomp; s/(.*)//; $name=$&; print ">$name$_" if tr/A-Z/A-Z/ >= 2000' assembly/spades/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/scaffolds.fasta > assembly/spades/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/scaffolds_gt2000.fasta

