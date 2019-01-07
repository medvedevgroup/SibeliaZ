import sys
import collections
MafRecord = collections.namedtuple('MafRecord', ['seq_name', 'start', 'size', 'strand', 'seq_size', 'body'])

def parse_maf(file_name):
	handle = open(file_name)
	alignment = []
	group = []
	for line in handle:
		line = line.strip()
		if line:
			if line[0] == 'a':
				alignment.append(group)
				group = []
			elif line[0] == 's':
				line = line.split()
				record = MafRecord(seq_name=line[1], start=int(line[2]), size=int(line[3]),
						 strand=line[4], seq_size=int(line[5]), body=line[6])
				group.append(record)
	alignment.append(group)
	alignment = [x for x in alignment if x]
	handle.close()
	return alignment


true_alignment = parse_maf(sys.argv[1])
min_size = min((len(x.body) for g in true_alignment for x in g if len(x.body)))
print min_size
#print [x for g in true_alignment for x in g if len(x.body) == min_size]

