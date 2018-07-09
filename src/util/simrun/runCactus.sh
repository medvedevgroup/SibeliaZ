#!/bin/bash
rm -rf ./result/cactus/$1
mkdir ./result/cactus/$1
python gen_cactus_input.py ./data/$1/ >> ./result/cactus/$1/cactus_input.txt

/research/ium125/progressiveCactus/bin/runProgressiveCactus.sh ./result/cactus/$1/cactus_input.txt \
./result/cactus/$1/tmp \
./result/cactus/$1/out.hal --maxThreads 64

source /research/ium125/progressiveCactus/environment && hal2mafMP.py ./result/cactus/$1/out.hal ./result/cactus/$1/cactus_unref.maf
python strip_cactus.py ./result/cactus/$1/cactus_unref.maf > ./result/cactus/$1/cactus.maf
./mafComparator --maf1=./data/$1/alignment.maf --maf2=./result/cactus/$1/cactus.maf --out=./result/cactus/$1/out.xml
