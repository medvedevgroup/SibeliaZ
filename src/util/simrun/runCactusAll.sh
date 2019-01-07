for i in 1 3 6 9 12 15 18 21 24
do
	/usr/bin/time -f "Cactus: %e elapsed %M memory" ./runCactus.sh $i 2> ./result/cactus/rt$i.txt
done
