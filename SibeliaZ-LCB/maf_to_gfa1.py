import sys
import numpy
import argparse
import collections
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
MafRecord = collections.namedtuple('MafRecord', ['seq_name', 'start', 'size', 'strand', 'seq_size', 'body'])

def pos_start(record):
	if record.strand == '+':
		return record.start
	return record.seq_size - (record.start + record.size)

def parse_maf(file_name):
        handle = open(file_name)
        group = []
        for line in handle:
                line = line.strip()
                if line:
                        if line[0] == 'a':
                                if len(group) > 0:
                                        yield group
                                group = []
                        elif line[0] == 's':
                                line = line.split()
                                record = MafRecord(seq_name=line[1], start=int(line[2]), size=int(line[3]),
                                                 strand=line[4], seq_size=int(line[5]), body=line[6])
                                group.append(record)
        handle.close()
        yield group

def print_maf(group, handle):
        for line in group:
                print >> handle, '\t'.join(str(x) for x in ('s', line.seq_name, line.start, line.size, line.strand, line.seq_size, line.body))

def print_all_maf(handle, alignment):
        print >> handle, '##maf version=1'
        for idx, group in enumerate(alignment):
                print >> handle, ''
                print >> handle, 'a'
                print_maf(group, handle)

def profile(group, column):
        return [group[i].body[column] == '-' for i in xrange(len(group))]

def slice(group, column_start, column_end, pos, shift):
        return [MafRecord(seq_name=group[i].seq_name, start=pos[i] + shift, size=column_end - column_start, strand=group[i].strand, seq_size=group[i].seq_size, body=group[i].body[column_start:column_end]) for i in xrange(len(group)) if group[i].body[column_start] != '-']

def is_homogeneous(group, column):
	ch = set()
	for record in group:
		ch.add(record.body[column])
	return len(ch) == 1 or (len(ch) == 2 and '-' in ch)

def decompose_column(group, blocks):
	group.sort(key=lambda x: x.body[0])
	i = 0
	while i < len(group):
		j = i
		while j < len(group) and group[i].body[0] == group[j].body[0]:
			j += 1
		blocks.append(group[i:j])
		i = j

# Split a MAF subrange with the same gap patter into subblocks with identical strings
def split_range(group, column_start, column_end, pos, blocks):
	origin = column_start
	for i in xrange(column_start, column_end):
		if not is_homogeneous(group, i):
			if column_start < i:
				blocks.append(slice(group, column_start, i, pos, column_start - origin))
			bad_column = slice(group, i, i + 1, pos, i - origin)
			decompose_column(bad_column, blocks)
			column_start = i + 1
	if column_start < column_end:
		blocks.append(slice(group, column_start, column_end, pos, column_start - origin))


# Split MAF blocks into subblocks where aligned strings are identical
# It works by first cutting MAF into subranges where the pattern of gaps is identical
# For example:
#	AC--T
#	AGTCT
#	ACTCT
# Will produces ranges (0, 1) (1, 2), (2, 4), (4, 5)
def split_maf_blocks(maf_file):
	blocks = []
	sequence = dict()
	for maf in parse_maf(maf_file):
		for record in maf:
			if record.seq_name not in sequence:
				sequence[record.seq_name] = []
	        pos = [record.start for record in maf]
        	prev_profile = profile(maf, 0)
	        prev_column = 0
        	while prev_column < len(maf[0].body):
                	next_column = prev_column
	                pos_inc = [0 for _ in xrange(len(pos))]
        	        while next_column < len(maf[0].body):
                	        next_profile = profile(maf, next_column)
                        	for i in xrange(len(pos)):
                                	pos_inc[i] += 0 if maf[i].body[next_column] == '-' else 1
	                        if next_profile == prev_profile:
        	                        next_column += 1
                	        else:
                        	        prev_profile = next_profile
                                	break

	                split_range(maf, prev_column, next_column, pos, blocks)
	                for i in xrange(len(pos)):
        	                pos[i] += pos_inc[i] - (0 if next_column < len(maf[0].body) and maf[i].body[next_column] == '-' else 1)
                	prev_column = next_column
	return (blocks, sequence)

# Generate blocks from the sequences uncovered by MAF
def get_uncovered_blocks(fasta, blocks, sequence):
	covered = dict()
	sequence_record = dict()
	for fasta_file in fasta:
	        for record in SeqIO.parse(fasta_file, "fasta"):
        	        sequence_record[record.id] = record
                	covered[record.id] = [False for _ in xrange(len(record.seq))]
			if record.id not in sequence:
				sequence[record.id] = []

        for b in xrange(len(blocks)):
               	for record in blocks[b]:
                       	sequence[record.seq_name].append((pos_start(record), b, record))
                        for i in xrange(pos_start(record), pos_start(record) + record.size):
       	                        covered[record.seq_name][i] = True

        for seq_id, cov in covered.items():
       	        i = 0
               	while i < len(cov):
                       	if cov[i] == False:
                               	j = i
                                while j < len(cov) and cov[j] == False:
       	                                j += 1
               	                blocks.append([MafRecord(seq_name=seq_id, start=i, size=j - i, strand='+', seq_size=len(cov), body=sequence_record[seq_id].seq[i:j])])
                       	        sequence[seq_id].append((i, len(blocks) - 1, blocks[-1][0]))
                               	i = j
                        else:
       	                        i += 1

def blocks_debug_output(blocks):
	maf_out = open("out.maf", "w")
	print_all_maf(maf_out, blocks)
	maf_out.close()

def output_block(b, remember_block):
        if b not in remember_block:
                print "S\t" + str(b + 1) + "\t" + blocks[b][0].body
                remember_block.add(b)

def output_link(a, b, remember_block, remember_link):
        start1, block1, record1 = a
        start2, block2, record2 = b
        output_block(block1, remember_block)
        output_block(block2, remember_block)
        if start1 + record1.size == start2:
                link = ','.join((str(block1), record1.strand, str(block2), record2.strand))
                if link not in remember_link:
                        id = len(remember_link)
                        remember_link[link] = id
                        print "\t".join(("L", str(block1 + 1), record1.strand, str(block2 + 1), record2.strand, "*"))
                id = remember_link[link]
        else:
                print "FAIL", start1, record1.size, start2

parser = argparse.ArgumentParser(description='A helper script for covnerting MAF produced by SibeliaZ to GFA1.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('maf', help='MAF output by SibeliaZ')
parser.add_argument('fasta', nargs='+', help='Input genomes')
args = parser.parse_args()

blocks, sequence = split_maf_blocks(args.maf)
get_uncovered_blocks(args.fasta, blocks, sequence)

#blocks_debug_output(blocks)

remember_block = set()
remember_link = dict()

print "H\tVN:Z:1.0"

for header, blocks_seq in sequence.items():
	blocks_seq.sort()
	for i in xrange(0, len(blocks_seq) - 1):
		output_link(blocks_seq[i], blocks_seq[i + 1], remember_block, remember_link)
	print "P\t" + header + "\t" + ','.join((str(block[1] + 1) + block[2].strand for block in blocks_seq))

