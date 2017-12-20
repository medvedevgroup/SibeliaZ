import sys
from Bio import SeqIO

class GenomeCollection:
	def __init__(self):
		self.seq = dict()
		self.order = []

	def Add(self, seq_species, seq_chr, seq_length):
		if seq_species not in self.seq:
			self.seq[seq_species] = dict()
		self.seq[seq_species][seq_chr] = len(self.order) + 1
		self.order.append((seq_species, seq_chr, seq_length))

	def WriteHeader(self):
		print('Seq_id\tSize\tDescription')
		for idx, (seq_species, seq_id, seq_length) in enumerate(self.order):
			print('\t'.join(str(x) for x in(idx, seq_species + "." + seq_id, seq_length)))
		print('--------------------------------------------------------------------------------')

	def WriteRecords(self, annotation):

class Annotation:
	def __init__(self):
		self.record=[]

	def AddHeader(self, header):
		self.header = dict()
		for val, key in enumerate(header):
			self.header[key] = val

	def AddRecord(self, rec):
		self.record.append(rec)

def parse_annotation(file_name):
	ret = Annotation()
	with open(file_name) as handle:
		for line in handle:
			if 'FBgn_ID' in line:
				ret.AddHeader(line.strip().split())
			if 'Dyak' in line:
				ret.AddRecord(line.strip().split())
	return ret

annotation = parse_annotation(sys.argv[2])
genomes = GenomeCollection()
with open(sys.argv[1], "rU") as handle:
	for record in SeqIO.parse(handle, "fasta"):
		seq_id = record.id
		seq_length = len(record.seq)
		description = record.description.split('.')
		seq_species, seq_chr = description[0], description[1]
		genomes.Add(seq_species, seq_chr, seq_length)

genomes.WriteHeader(annotation)