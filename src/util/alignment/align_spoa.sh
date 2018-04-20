#find ./gene_blocks/ -name "*.fa" -printf "%f\n" 
find ./gene_blocks/ -name "*.fa" -printf "%f\n" | xargs -I @ -P 63 ./run_spoa.sh @
