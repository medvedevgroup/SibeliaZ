#!/bin/bash
rm -rf ./result/lsibelia/$1
mkdir ./result/lsibelia/$1
k=$2
tr=1
rt=./result/lsibelia/rt$1.txt

/usr/bin/time -f "%e elapsed\n%M memory" ./twopaco -k $k -f 32 -t 1 -o ./result/lsibelia/$1/$1_$k.dbg ./data/$1/genomes.fa 2> $rt

/usr/bin/time -f "%e elapsed\n%M memory" ./sibelia-pp \
--gfile ./data/$1/genomes.fa \
--infile ./result/lsibelia/$1/$1_$k.dbg \
-o ./result/lsibelia/$1/$1_$k \
-k $k -b 150 -m 100 -t $tr --abundance 512 2>> $rt

rm -rf ./align
mkdir ./align
find ./result/lsibelia/$1/$1_$k/blocks -name "*.fa" -printf "%p\n" | xargs -I @ -P 8 ./run_mlagan.sh @
python collect_lagan.py ./align > ./result/lsibelia/$1/$k.maf
./mafComparator --maf1=./data/$1/alignment.maf --maf2=./result/lsibelia/$1/$k.maf --out=./result/lsibelia/$1/out.xml
rm -rf ./align
