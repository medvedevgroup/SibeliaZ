import sys
import math
import collections
from Bio import SeqIO
from os import listdir
from os.path import isfile, join
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

gene = []
gene_id = dict()

mypath = './gene_seq/'
for f in listdir(mypath):
	p = join(mypath, f)
	if isfile(p):
	        for record in SeqIO.parse(p, "fasta"):
			gene.append(record)
			gene_id[record.id] = len(gene) - 1

block_n = 0
mypath = './homolog'
for f in listdir(mypath):
        p = join(mypath, f)
	seen = set()
        if isfile(p) and ".txt" in f:
		for line in open(p):
			line = line.strip().split()
			if len(line) < 2:
				continue
			id = [gene_id[x] for x in line[:2] if x in gene_id]
			if len(id) == 2:
				check = line[:2]
				check.sort()
				check = '\t'.join(check)
				if check in seen:
					print check
					continue
				seen.add(check)
				score = round(max(float(line[2]), float(line[3])))
				file = open("./gene_blocks/" + str(block_n) + "_" + str(int(score)) + ".fa", "w")
				for gid in id:
					SeqIO.write(gene[gid], file, "fasta")
				block_n += 1
