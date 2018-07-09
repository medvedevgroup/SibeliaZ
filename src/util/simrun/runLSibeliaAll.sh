rm -rf ./result/lsibelia
mkdir ./result/lsibelia

for i in 1 3 6 9 12 15 18 21 24
do
	rtdir=/research/ium125/simrun/result/lsibelia/$i
	mkdir $rtdir
	for k in {15..35..2}
	do
		kdir=$rtdir/$k
		mkdir $kdir
		for b in {0..200..25}
		do
			for f in {0..50..25}
			do
				nwdir=$kdir/$b\_$f
				mkdir $nwdir
				./runLSibelia.sh $i $k $b $f $nwdir 
			done
		done
	done
done
