import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

for record in SeqIO.parse(sys.argv[1], "fasta"):
	new_seq = ''.join(("N" if x.islower() else x for x in record.seq))
	seq_species = record.description.split()[-1].split('=')[-1][:-1]
	new_record = SeqRecord(Seq(new_seq, IUPAC.ambiguous_dna), id=seq_species + "." + record.id, description="")
	SeqIO.write(new_record, sys.stdout, "fasta")
	
	