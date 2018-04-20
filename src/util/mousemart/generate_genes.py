import sys
import collections
from os import listdir
from os.path import isfile, join
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord


chr_name_to_genbank = dict()
missing = []

mypath = './assembly_reports/'
for f in listdir(mypath):
	p = join(mypath, f)
	if isfile(p):
		strain = f.split('.')[1]
		strain = 'm' + ''.join(strain.split('_')[1:-3]).lower()
		chr_name_to_genbank[strain] = dict()
		for line in open(p):
			if line[0] != '#':
				line = line.strip().split()
				chr_name_to_genbank[strain][line[0]] = line[4]

strain_chr_seq = dict()
mypath = './genome/'
for f in listdir(mypath):
        p = join(mypath, f)
        if isfile(p):
		strain = f.split('.')[1]
                strain = 'm' + ''.join(strain.split('_')[1:-2]).lower()
		strain_chr_seq[strain] = dict()
		for record in SeqIO.parse(p, "fasta"):
    			strain_chr_seq[strain][record.id] = record.seq

mypath = './gene/'
print >> sys.stderr, chr_name_to_genbank
for f in listdir(mypath):
        p = join(mypath, f)
	if '.txt' in f:
		strain = f.split('_')[0]
		if not strain in strain_chr_seq:
			continue
		for line in open(p):
			line = line.strip().split('\t')
			id, chr, start, end, strand = line
			start = int(start) - 1
			end = int(end)
			seq_length = end - start
			if chr not in chr_name_to_genbank[strain]:
				print >> sys.stderr, id, chr
				continue
			seq_genbank_id = chr_name_to_genbank[strain][chr]
			seq = strain_chr_seq[strain][seq_genbank_id][start:end]
			length = len(strain_chr_seq[strain][seq_genbank_id])
			if strand == '-1':
				seq = seq.reverse_complement()
				start = length - end
			descr = ';'.join((seq_genbank_id, str(start), str(seq_length), strand, str(length)))
			rec = SeqRecord(Seq(str(seq), generic_dna), id=id, description=descr)
			file = open("./gene_seq/" + id + ".fa", "w")
			SeqIO.write(rec, file, 'fasta')
			file.close()

