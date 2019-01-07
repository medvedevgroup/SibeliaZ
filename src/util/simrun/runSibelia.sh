#!/bin/bash

k=$2
b=$3
dir=$4
rt=$dir/rt$1.txt

echo "1" > run.txt
echo "$k $b" >> run.txt

/usr/bin/time -f "Sibelia: %e elapsed %M memory" ./Sibelia \
./data/$1/genomes.fa \
-o $dir \
-k run.txt --lastk 100 -m 100 --nopostprocess 2> $rt

/usr/bin/time -f "Alignment: %e elapsed %M memory" python C-Sibelia.py -o $dir ./data/$1/genomes.fa -p 16 -a $dir/$k.maf 2>> $rt
./mafComparator --maf1=./data/$1/alignment.maf --maf2=$dir/$k.maf --out=$dir/out.xml

