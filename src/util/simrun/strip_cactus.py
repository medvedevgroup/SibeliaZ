import sys

for line in open(sys.argv[1]):
	orig = line
	line = line.strip().split('\t')
	if line[0] == 's' and 'SE' in orig:
		line[1] = line[1].split('.')[0]
	print '\t'.join(line)
