import sys
import collections
from Bio import SeqIO
from os import listdir
from os.path import isfile, join
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import connected_components

gene = []
gene_id = dict()

mypath = './gene_seq/'
for f in listdir(mypath):
        p = join(mypath, f)
#	if len(gene) > 100:
#		break
        if isfile(p):
	        for record in SeqIO.parse(p, "fasta"):
			gene.append(record)
			gene_id[record.id] = len(gene) - 1

n = len(gene)
graph = lil_matrix((n, n))

block_n = 0
mypath = './homolog'
for f in listdir(mypath):
        p = join(mypath, f)
        if isfile(p) and ".txt" in f:
		for line in open(p):
			id = [gene_id[x] for x in line.strip().split() if x in gene_id]
			if len(id) > 1:
				file = open("./gene_blocks/" + str(block_n) + ".fa", "w")
				for gid in id:
					SeqIO.write(gene[gid], file, "fasta")
				block_n += 1
