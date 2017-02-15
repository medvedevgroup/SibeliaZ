import sys
import collections
from Bio import SeqIO

for record in SeqIO.parse(sys.argv[1], "fasta"):
	print '>' + record.id
	print ''.join((x for x in record.seq if x != '-'))


