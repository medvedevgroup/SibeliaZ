import sys
from os import listdir
from os.path import isfile, join


missing = 0
mypath = './alignment/'
handle = [open('alignment_bin_chr_y/' + str(i) + '.maf', 'w') for i in xrange(0, 10)]
for h in handle:
	print >> h, "##maf version=1"


for f in listdir(mypath):
        p = join(mypath, f)
        if isfile(p):
		h = open(p)
		line = [line.strip() for line in h]
		if len(line) > 1:
			al1 = line[1].split()[-1]
			al2 = line[2].split()[-1]
			if ("CM001014.2" in line[1]) and ("CM001014.2" in line[2]):
				match = len([i for i in xrange(len(al1)) if al1[i] == al2[i]])
				identity1 = float(match) / int(line[1].split()[3])
				identity2 = float(match) / int(line[2].split()[3])
				identity = max(identity1, identity2)
				bin = min(9, int(identity * 10))
				print >> handle[bin], ''
	                        print >> handle[bin], '\n'.join(line)
		else:
			missing.append(p)
		h.close()

print >> sys.stderr, missing

