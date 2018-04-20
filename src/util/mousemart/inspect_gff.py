import sys
import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

MafRecord = collections.namedtuple('MafRecord', ['seq_name', 'strand', 'start', 'size', 'seq_size', 'body'])

def parse_maf(file_name):
	group = []
	handle = open(file_name)
	for line in handle:
		line = line.strip()
		if line:
			if line[0] == 'a':
				if group:
					yield group
					group = []
			elif line[0] == 's':
				line = line.split()
				record = MafRecord(seq_name=line[1], start=int(line[2]), size=int(line[3]),
						 strand=line[4], seq_size=int(line[5]), body=line[6])
				group.append(record)
	if group:
		yield group

count = 0
covered = dict()
for record in SeqIO.parse(sys.argv[1], "fasta"):
	covered[record.id] = [False] * len(record.seq)

for line in open(sys.argv[2]):
	if not line or line[0] == '#':
		continue
	line = line.strip().split()
	seq_name = line[0]
	start, end = int(line[3]), int(line[4])
	covered[seq_name][start:end] = [True] * (end - start)

def print_maf(group):
	print 'a'
	for line in group:
		print '\t'.join((str(x) for x in ('s', line.seq_name, line.start, line.size, line.strand, line.seq_size, line.body)))

for record in parse_maf(sys.argv[3]):
	total_length = 0
	covered_length = 0
        for line in record:
		total_length += line.size
                if line.seq_name not in covered:
			print "Unknown sequence:", line.seq_name
			continue
                start = line.start
                if line.strand == '-':
                        start = line.seq_size - (start + line.size)
		covered_length += len([x for x in covered[line.seq_name][start:start + line.size] if x])
	ratio = float(covered_length) / total_length
	if ratio < 0.8:
		print ratio
		print_maf(record)

