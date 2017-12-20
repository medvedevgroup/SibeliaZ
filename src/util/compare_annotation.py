import sys
from Bio import SeqIO

class GenomeCollection:
	def __init__(self):
		self.seq = dict()
		self.order = []

	def Add(self, seq_species, seq_chr, seq_length):
		if seq_species not in self.seq:
			self.seq[seq_species] = dict()
		self.seq[seq_species][seq_chr] = str(len(self.order) + 1)
		self.order.append((seq_species, seq_chr, seq_length))

	def SeqId(self, species, chrom):
		return self.seq[species][chrom]

	def WriteHeader(self):
		print('Seq_id\tSize\tDescription')
		for idx, (seq_species, seq_id, seq_length) in enumerate(self.order):
			print('\t'.join(str(x) for x in(idx + 1, seq_length, seq_species + "." + seq_id)))
		print('--------------------------------------------------------------------------------')

class Annotation:
	def __init__(self):
		self.record=[]

	def AddHeader(self, header):
		self.header = dict()
		for val, key in enumerate(header):
			self.header[key] = val - 1

	def AddRecord(self, rec):
		self.record.append(rec)

	def WriteBlock(self, genomes, species, chrom, strand, start, end):
		length = str(int(end) - int(start))
		if strand > 0:
			print('\t'.join((genomes.SeqId(species, chrom), "+", start, end, length)))
		else:
			print('\t'.join((genomes.SeqId(species, chrom), "-", end, start, length)))		

	def WriteRecords(self, genomes):
		src_species = 'Dmel'
		for block_cnt, record in enumerate(self.record):
			print("Block #" + str(block_cnt))
			print("Seq_id\tStrand\tStart\tEnd\tLength")
			src_chr = record[self.header['Arm/Scaffold']]			
			src_strand = int(record[self.header['Strand']])
			src_start, src_end =  record[self.header['Location']].split('..')
			self.WriteBlock(genomes, src_species, src_chr, src_strand, src_start, src_end)
			target_species = record[self.header['Ortholog_GeneSymbol']].split('\\')[0]
			target_chr = record[self.header['Ortholog_Arm/Scaffold']]
			target_strand = int(record[self.header['Ortholog_Strand']])
			target_start, target_end = record[self.header['Ortholog_Location']].split('..')
			self.WriteBlock(genomes, target_species, target_chr, target_strand, target_start, target_end)
			print('--------------------------------------------------------------------------------')

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

genomes.WriteHeader()
annotation.WriteRecords(genomes)