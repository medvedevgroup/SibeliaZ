import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

for record in SeqIO.parse(sys.argv[1], "fasta"):
	new_seq = ''.join(("N" if x.islower() else x for x in record.seq))
	new_record = SeqRecord(Seq(new_seq, IUPAC.ambiguous_dna), id=record.id, name=record.name, description=record.description)
	SeqIO.write(new_record, sys.stdout, "fasta")
	
	