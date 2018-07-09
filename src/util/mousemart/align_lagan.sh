find ./gene_blocks/ -name "*.fa" -printf "%f\n" | xargs -I @ -P 16 ./run_lagan.sh @
