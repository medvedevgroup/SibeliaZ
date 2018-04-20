import sys
import collections
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
for record in parse_maf(sys.argv[1]):
	count += 1
	if count % 100 == 0:
		print >> sys.stderr, count
	for line in record:
		if line.seq_name not in covered:
			covered[line.seq_name] = [False for _ in xrange(line.seq_size)]
		start = line.start
		if line.strand == '-':
			start = line.seq_size - (start + line.size)
		covered[line.seq_name][start:start + line.size] = [True] * line.size

def print_maf(group):
	print 'a'
	for line in group:
		print '\t'.join((str(x) for x in ('s', line.seq_name, line.start, line.size, line.strand, line.seq_size, line.body)))

for record in parse_maf(sys.argv[2]):
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

