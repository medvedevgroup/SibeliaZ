#!/bin/bash
k=$2
b=$3
f=$4
dir=$5
tr=1
rt=$dir/rt$1.txt

/usr/bin/time -f "TwoPaco: %e elapsed %M memory" ./twopaco -k $k -f 32 -t 1 -o $dir/graph.dbg ./data/$1/genomes.fa 2> $rt

/usr/bin/time -f "L-Sibelia: %e elapsed %M memory" ./sibelia-pp \
--gfile ./data/$1/genomes.fa \
--infile $dir/graph.dbg \
-o $dir \
-k $k -b $b -m 100 -t $tr -f $f 2>> $rt

/usr/bin/time -f "Alignment: %e elapsed %M memory" python C-Sibelia.py -o $dir ./data/$1/genomes.fa -p 16 -a $dir/$k.maf 2>> $rt
./mafComparator --maf1=./data/$1/alignment.maf --maf2=$dir/$k.maf --out=$dir/out.xml
