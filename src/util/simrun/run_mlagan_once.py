#!/usr/local/bin/python2.7

import os
import sys
import time
import shutil
import tempfile
import argparse
import itertools
import functools
import subprocess
import collections
import multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord


INSTALL_DIR = os.path.dirname(os.path.abspath(__file__))
LAGAN_DIR = os.path.join(INSTALL_DIR, '..', 'lib', 'Sibelia', 'lagan')
os.environ['LAGAN_DIR'] = LAGAN_DIR

SyntenyBlock = collections.namedtuple('SyntenyBlock', ['seq', 'chr_id', 'strand', 'id', 'start', 'end', 'chr_num', 'chr_size'])
AlignmentRecord = collections.namedtuple('AlignmentRecord', ['body', 'block_instance'])


def process_block(block_file, outfile):
	tmp_dir = tempfile.mkdtemp()
	instance_list = [record for record in SeqIO.parse(block_file, "fasta")]
	os.chdir(tmp_dir)
	alignment_file = "al.fa"
	file_name = [tmp_dir + "/" + str(i) for i, _ in enumerate(instance_list)]
	cmd = [os.path.join(LAGAN_DIR, "mlagan")] + file_name
	alignment_handle = open(alignment_file, 'w')
	for index, block in enumerate(instance_list):
		SeqIO.write(block, file_name[index], "fasta")
	worker = subprocess.Popen(cmd, stdout=alignment_handle, stderr=subprocess.PIPE)
	_, stderr = worker.communicate()
	if worker.returncode != 0:
		print >> sys.stderr, "Error on block:", block_file, stderr
	alignment_handle.close()
	print >> sys.stdout, "a"
	for record in SeqIO.parse(alignment_file, "fasta"):
		description = record.description.split(";")
		description[0] = description[0].split()[1]
                description[-2] = '+' if description[-2] == '1' else '-'
                data = ["s"] + description + [str(record.seq)]
                print >> sys.stdout, '\t'.join(data)
	shutil.rmtree(tmp_dir)


start = time.time()
parser = argparse.ArgumentParser(description='A tool for comparing two microbial genomes.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('block', help='Output of L-Sibelia')
args = parser.parse_args()
process_block(args.block, "")
