import sys
from Bio import SeqIO

class GenomeCollection:
	def __init__(self):
		self.seq = dict()
	def Add(self, seq_id, seq_length, seq_species):
		if seq_species not in self.seq:
			self.seq[seq_species] = dict()
		self.seq[seq_species][seq_id] = [0 for _ in xrange(seq_length)]

genomes = GenomeCollection()
with open(sys.argv[1], "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
	seq_id = record.id
	seq_length = len(record.seq)
	seq_species = record.description.split()[-1].split('=')[-1][:-1]
	genomes.Add(seq_id, seq_length, seq_species)