#!/usr/local/bin/python2.7

import re
import os
import sys
import time
import glob
import shutil
import tempfile
import argparse
import itertools
import functools
import subprocess
import collections
import multiprocessing

COVER = 1
UNCOVER = 0
LINE_LENGTH = 60
INSTALL_DIR = os.path.dirname(os.path.abspath(__file__))

FastaRecord = collections.namedtuple('FastaRecord', ['seq', 'description', 'id'])
SyntenyBlock = collections.namedtuple('SyntenyBlock', ['seq', 'chr_id', 'strand', 'id', 'start', 'end', 'chr_num', 'chr_size'])
AlignmentRecord = collections.namedtuple('AlignmentRecord', ['body', 'block_instance'])

class FailedStartException(Exception):
	pass

class DuplicatedSequenceIdException(Exception):
	def __init__(self, id):
		self._id = id

	def __str__(self):
		return 'Found duplicated sequence id "%s"' % (self._id)

def unzip_list(zipped_list):
	return ([x for (x, _) in zipped_list], [y for (_, y) in zipped_list])

def unzip_list3(zipped_list):
	return ([x for (x, _, _) in zipped_list], [y for (_, y, _) in zipped_list], [z for (_, _, z) in zipped_list])

def parse_fasta_file(file_name):
	handle = open(file_name)
	line = [line.strip() for line in handle if line.strip() != '']
	record = []
	i = 0
	while i < len(line):
		if line[i][0] == '>':
			j = i + 1
			while j < len(line) and line[j][0] != '>':
				j += 1
			seq = ''.join(line[i + 1:j])
			description = line[i][1:].strip()
			seq_id = description.split()[0]
			record.append(FastaRecord(seq=seq, description=description, id=seq_id))
			i = j
		else:
			i += 1
	handle.close()
	return record

def reverse_complementary(seq):
	comp = dict()
	comp['A'] = 'T'
	comp['T'] = 'A'
	comp['C'] = 'G'
	comp['G'] = 'C'
	return ''.join([(comp[ch] if ch in comp else ch) for ch in seq[::-1]])

def parse_blocks_coords(blocks_file, genome):
	group = [[]]
	num_seq_size = dict()
	num_seq_id = dict()
	seq_id_num = dict()
	line = [l.strip() for l in open(blocks_file) if l.strip()]
	for l in line:
		if l[0] == '-':
			group.append([])
		else:
			group[-1].append(l)
	for l in group[0][1:]:
		l = l.split()
		num_seq_id[l[0]] = l[2]
		num_seq_size[int(l[0])] = int(l[1])
		seq_id_num[l[2]] = int(l[0])
	ret = dict()
	for g in [g for g in group[1:] if g]:
		block_id = int(g[0].split()[1][1:])
		ret[block_id] = []
		for l in g[2:]:
			l = l.split()
			chr_id = num_seq_id[l[0]]
			start = int(l[2])
			end = int(l[3])
			chr_num = int(l[0])
			strand = l[1]
			if strand == '+':
				true_start = start - 1
				true_end = end
			else:
				true_start = end - 1
				true_end = start
			seq_num = int(l[0]) - 1
			seq = genome[seq_num].seq[true_start:true_end]
			if strand == '-':
				seq = reverse_complementary(seq)
			ret[block_id].append(SyntenyBlock(seq=seq, chr_id=chr_id, strand=strand, id=block_id, start=start,
							end=end, chr_num=chr_num, chr_size=num_seq_size[chr_num]))
	return ret

def strip_chr_id(chr_id):
	part = chr_id.split('|')
	if len(part) == 5:
		return part[-2].split('.')[0]
	return chr_id

def write_wrapped_text(text, handle):
	pos = 0
	while pos < len(text):
		end = min(pos + LINE_LENGTH, len(text))
		print >> handle, text[pos:end]
		pos = end

def write_fasta_records(fasta_record, file_name):
	handle = open(file_name, 'w')
	for record in fasta_record:
		print >> handle, '>' + record.description
		write_wrapped_text(record.seq, handle)
	handle.close()


def no_gaps(sequence):
	return ''.join([ch for ch in sequence if ch != '-'])

def get_seq(file_name):
	all_seq = [record for record in parse_fasta_file(file_name)]
	return [(record.id, record.seq) for record in all_seq ]

def parse_header(header):
	ret = dict()
	header = header.split(',')
	for item in header:
		item = item.split('=')
		key = item[0]
		value = ''.join([ch for ch in item[1] if ch != "'" and ch != '"'])
		ret[key] = value
	return ret

def find_instance(instance_list, reference_seq_id, in_reference):
	for instance in instance_list:
		if (instance.chr_id in reference_seq_id) == in_reference:
			return instance
	return None

def process_block(block):
	pid = str(os.getpid()) + '_'
	alignment_file = pid + 'align.fasta'
	synteny_block_id, instance_list = block
	file_name = pid + 'block.fasta'
	spoa_cmd = ['../spoa', file_name, '-l', '1', '-r', '1']
	fasta_records = []
	for index, block in enumerate(instance_list):
		fasta_records.append(FastaRecord(id=block.chr_id, description=str(index), seq=block.seq))

	write_fasta_records(fasta_records, file_name)
	worker = subprocess.Popen(spoa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = worker.communicate()
	if worker.returncode != 0:
		raise FailedStartException(stderr)

	alignment = stdout.strip().split('\n')[1:]
	ret = [AlignmentRecord(body=align, block_instance=instance_list[idx]) for idx, align in enumerate(alignment)]
	os.remove(file_name)
	return ret

def call_variants(directory, genomes, proc_num):
	os.chdir(directory)
	block_seq = parse_blocks_coords('blocks_coords.txt', genomes)
	pool = multiprocessing.Pool(proc_num)
	annotated_block = []
	for synteny_block_id, instance_list in block_seq.items():
		annotated_block.append((synteny_block_id, instance_list))
	if annotated_block:
		alignment = pool.map_async(process_block, annotated_block).get()
		pool.close()
		pool.join()
	else:
		alignment = []

	os.chdir('..')
	return alignment


def write_alignments_maf(alignment_list, noparalog, handle):
	print >> handle, '##maf version=1\n'
	for group in alignment_list:
		print >> handle, 'a'
		for alignment in group:
			block = alignment.block_instance
			start = min(block.start, block.end) - 1
			end = max(block.start, block.end)
			if block.strand != '+':
				start = block.chr_size - end
			print >> handle, 's', block.chr_id, start, abs(block.end - block.start) + 1, block.strand, block.chr_size, alignment.body
		print >> handle, ''

def handle_exception(e, temp_dir):
	print 'An error occured:', e
	shutil.rmtree(temp_dir, ignore_errors=True)

start = time.time()
parser = argparse.ArgumentParser(description='Batch aligner for Sibelia++', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('genome', nargs='*', help='FASTA files')
parser.add_argument('-p', '--processcount', help='Number of running processes', type=int, default=1)
parser.add_argument('--debug', help='Generate output in text files', action='store_true')
group = parser.add_mutually_exclusive_group()
group.add_argument('-o', '--outdir', help='Directory for synteny block output files')
args = parser.parse_args()

try:
	if args.outdir is None:
		if args.tempdir is None:
			try:
				temp_dir = tempfile.mkdtemp(dir='.')
			except EnvironmentError as e:
				print e
				exit(1)
		else:
			temp_dir = args.tempdir
	else:
		temp_dir = args.outdir

	genomes = list(itertools.chain.from_iterable([parse_fasta_file(file_name) for file_name in args.genome]))


	print >> sys.stderr, "Performing alignment..."
	alignment_list = call_variants(temp_dir, genomes, args.processcount)
	handle = open('alignment.maf', 'w')
	write_alignments_maf(alignment_list, False, handle)
	handle.close()

except FailedStartException as e:
	handle_exception(e, temp_dir)
except EnvironmentError as e:
	handle_exception(e, temp_dir)
except DuplicatedSequenceIdException as e:
	handle_exception(e, temp_dir)



