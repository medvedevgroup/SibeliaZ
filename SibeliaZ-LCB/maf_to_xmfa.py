import sys
from Bio import AlignIO

input_handle = sys.stdin
output_handle = sys.stdout

alignments = AlignIO.parse(input_handle, "maf")
AlignIO.write(alignments, output_handle, "mauve")

output_handle.close()
input_handle.close()
