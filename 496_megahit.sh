mkdir -p assembly/megahit
~/apps/megahit/megahit --12 reads/normalized/ASGARD_20181114/all_pairs.fastq.gz -r reads/normalized/ASGARD_20181114/all_singles.fastq.gz -r reads/normalized/ASGARD_20181114/all_merged.fastq.gz -t 7 -m 512e9 -o assembly/megahit/all_ASGARD_20181114_bio --out-prefix ASGARD_20181114 --min-contig-len 1
#~/apps/megahit/megahit -o assembly/megahit/all_ASGARD_20181114 --out-prefix ASGARD_20181114 --continue

