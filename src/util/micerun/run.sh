echoerr() { echo "$@" 1>&2; }
for d in 2 4 8 16
do
	gdir="/research/ium125/Mousemart_new/genome"
	rdir="/research/ium125/TwoPaCo/build/graphconstructor"
	gfile=$gdir/gmice_$d.fa
	infile=gmice_$d\_29.dbg
	outdir=out_$d
	rm -rf outdir
	k=29
	echoerr "-------------------------------------------------"
	/usr/bin/time -f "TwoPaco: %e elapsed %M memory" ./twopaco -t 16 -k $k -f 34 $gfile -o $infile
	for a in 25 50 100 200 250
	do
		/usr/bin/time -f "L-Sibelia: %e elapsed %M memory" ./sibelia-pp --gfile $gfile --infile $infile -k 29 -b 200 -o $outdir\_$a -m 250 -t 16 --abundance $a
	done
done

