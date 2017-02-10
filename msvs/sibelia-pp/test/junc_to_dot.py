import sys

if len(sys.argv) != 2:
	print "Usage: junc_to_dot.py <junc_file>"
	exit()

prev_seq_id = -1
prev_junc_id = -1

print "digraph G\n{\nrankdir=LR"

for line in open(sys.argv[1]):
	seq_id, pos, junc_id = (int(x) for x in line.strip().split())
	if seq_id == prev_seq_id:
		print prev_junc_id, '->', junc_id, '[label="' + str(seq_id) + " " + str(pos) + '"]'
	prev_seq_id = seq_id
	prev_junc_id = junc_id

print "}"