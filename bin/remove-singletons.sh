#/bin/bash
file=$1

perl -0076 -ne '{chomp; m/(.*):(.*);/; print ">$_" if $1 ne $2};' $file

