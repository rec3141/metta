#for dir in assembly/spades/*; do if [ -d "$dir/K21" ]; then echo "keep $dir";  else echo "kill $dir"; rm -rf "$dir"; fi; done
#for dir in assembly/spades/*; do if [ -e $dir/corrected/*R1_001.fastq.00.0_0.cor.fastq.gz ]; then echo "keep $dir";  else echo "kill $dir"; ls -l $dir/corrected/*R1_001.fastq.00.0_0.cor.fastq.gz; echo rm -rf "$dir"; fi; done
for dir in assembly/spades/*; do 
	if [ -e $dir/corrected/*R1_001.fastq.00.0_0.cor.fastq.gz ]; then
#		echo "keep $dir";
#		echo "";
continue
	else 
		echo "kill $dir"; 
		ls -l $dir/corrected; 
		echo rm -rf "$dir"; 
	fi; 
done

