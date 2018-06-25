import sys
from os import listdir
from os.path import isfile, join


missing = 0
mypath = './alignment/'
handle = [open('alignment_bin_all/' + str(i) + '.maf', 'w') for i in xrange(0, 101)]

for f in listdir(mypath):
        p = join(mypath, f)
        if isfile(p):
		h = open(p)
		line = [line.strip() for line in h]
		if len(line) > 1:
			pam = int(f.split('.')[0].split('_')[1])
			print >> handle[pam], ''
			print >> handle[pam], '\n'.join(line)
		else:
			missing.append(p)
		h.close()

print >> sys.stderr, missing

