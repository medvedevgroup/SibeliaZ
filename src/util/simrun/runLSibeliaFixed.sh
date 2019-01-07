rm -rf ./result/lsibelia
mkdir ./result/lsibelia

for i in 1 3 6 9 12 15 18 21 24
do
	rtdir=/research/ium125/simrun/result/lsibelia/$i
	mkdir $rtdir
	k=15
	kdir=$rtdir/$k
	mkdir $kdir
	b=200
	f=25
	nwdir=$kdir/$b\_$f
	mkdir $nwdir
	./runLSibelia.sh $i $k $b $f $nwdir 
done
