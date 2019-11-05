#/bin/bash

module load lang/Anaconda3/2.5.0
source activate `pwd`/prokka

export PERL5LIB=/home/recollins/work/nanobase/Arctic_metagenomes/prokka/lib/site_perl/5.26.2:/home/recollins/work/nanobase/Arctic_metagenomes/prokka/lib/5.26.2

prokka --help
