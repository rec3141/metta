awk -F' ' '{if ($0 ~ />/) {print $1 "_" (NR+1)/2} else {print $0}}' plass-assembly.faa > plass-assembly2.faa &
