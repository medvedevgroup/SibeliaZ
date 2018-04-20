import sys
import collections
from os import listdir
from os.path import isfile, join
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord


strain_chr_seq = dict()
mypath = './genome/'

used_by_sibelia = dict()

for f in listdir(mypath):
        p = join(mypath, f)
        if isfile(p):
                for record in SeqIO.parse(p, "fasta"):
                	used_by_sibelia[record.id] = [False for _ in xrange(len(record.seq))]


