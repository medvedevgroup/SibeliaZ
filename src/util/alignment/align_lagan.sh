find /research/ium125/Sibelia-pp/builc/gmice_4_27_200_n/blocks -name "*.fa" -printf "%p\n" | xargs -I @ -P 61 ./run_mlagan.sh @

