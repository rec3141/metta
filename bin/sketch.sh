for file in assembly/spades/spades_Plate*/scaffolds.fasta; do head -c1000000 $file > out/tmpseq.fa; ~/apps/bbmap/sendsketch.sh in=out/tmpseq.fa address=nt; done

