rm -rf ./result/sibelia
mkdir ./result/sibelia

for i in 1 3 6 9 12 15 18 21 24
do
	rtdir=/research/ium125/simrun/result/sibelia/$i
	mkdir $rtdir
	for k in {15..35..2}
	do
		kdir=$rtdir/$k
		mkdir $kdir
		for b in {0..200..25}
		do
			nwdir=$kdir/$b
			mkdir $nwdir
			./runSibelia.sh $i $k $b $nwdir 
		done
	done
done
