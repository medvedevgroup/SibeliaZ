for d in alignment alignment_bin_all gene_seq gene gene_blocks homolog
do
	rm -rf $d
	mkdir $d
done

python gene_xml.py
./get_all.sh
python generate_genes.py
python generate_blocks_pwise.py
./align_lagan.sh
python collect_lagan_bin.py
