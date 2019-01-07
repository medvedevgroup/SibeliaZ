#!/bin/bash

for f in ./gene/*.xml
do
	./get.sh $f "$f.txt"
done

for f in ./homolog/*.xml
do
        ./get.sh $f "$f.txt"
done
