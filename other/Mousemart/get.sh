#!/bin/bash

q=`<$1`
wget -O $2 "http://www.ensembl.org/biomart/martservice?query=$q"
