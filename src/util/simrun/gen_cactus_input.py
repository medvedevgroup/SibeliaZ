import os
import sys

directory = sys.argv[1]
tree_path = os.path.join(directory, "MuRealTree.nwk")
if os.path.exists(tree_path):
	tree = open(tree_path).next().strip().split(':')
	print ':'.join(tree[:-1]) + ';'

for filename in os.listdir(directory):
	if filename.endswith(".fasta"):
		path = os.path.join(directory, filename)
		print(filename[:-6] + " " + os.path.abspath(path))
