#/bin/bash
#input $1 is directory DIR of fastq.gz files to rename
# this takes the part of the name up to the _L001 and changes _ to -
# outputs links to DIR_rename

DIR=$1

for file in $DIR/*.fastq.gz; do
	echo $file
	D1=`basename $DIR`
	F1=`basename $file`

	head=`echo $F1 | awk 'match($0,/L001/) {print substr($0,0,RSTART-2)}' | tr '_' '-'`
	tail=`echo $F1 | awk 'match($0,/L001/) {print substr($0,RSTART-1,999)}'`
	echo $head$tail
	ln -s `pwd`/$file `pwd`/$DIR"_rename"/$head$tail

done
