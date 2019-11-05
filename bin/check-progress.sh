#/bin/bash
#input $1 is directory to check
#for dir in `ls -d */ | cut -f1 -d'/' | sort -r`; do ./check-progress.sh $dir; done
for file in $1/tsv/*.tsv; do
	ls -l $file
	profile="`tail -n12 $file | head -n1 | tr -s ' ' | cut -f4 -d' '`$"
	grep ACC $1.hmm | awk -v var=$profile 'BEGIN{line=0}{if($0 ~ var) {line=NR}}; END{print 100*line/NR,line,NR};'
	tail $file | grep 'ok'
done
