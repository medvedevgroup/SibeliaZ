#!/bin/bash
find ./alignment_bin_chr_y/ -name "*.maf" -printf "%f\n" | xargs -I @ -P 63 ./compare_once.sh @
#find ./alignment_bin/ -name "*.maf" -printf "%f\n" 
