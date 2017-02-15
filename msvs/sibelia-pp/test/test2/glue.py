import sys
import collections
from Bio import SeqIO

record = []

for file_name in sys.argv[1:]:
           record.append([r for r in SeqIO.parse(file_name, "fasta")])

minlen = min(len(r) for r in record)

for i in xrange(0, minlen):
	print '>S' + str(i)
	print ''.join((str(r[i].seq) for r in record))


