#!/bin/bash

q=`<$1`
wget -O result.txt "http://www.ensembl.org/biomart/martservice?query=$q"
