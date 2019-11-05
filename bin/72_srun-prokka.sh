#!/bin/bash
#input $1 is assembly fasta file
#input $2 is metabat depth file

assembly=$1
metabat=$2

sbatch --partition=bio --ntasks=7 --tasks-per-node=7 --mem=1400G <<-EOF
#!/bin/bash
bin/71_run-prokka.sh $assembly $metabat
EOF

