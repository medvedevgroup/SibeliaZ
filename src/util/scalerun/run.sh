for t in 20 16 12 8 4 2 1
do
	/usr/bin/time -f "%e elapsed\n%M memory" ./sibelia-pp --gfile gmice_2.fa --infile gmice_2_29.dbg -k 29 -f -25 -b 200 -o out -m 250 -t $t --abundance 128
done

