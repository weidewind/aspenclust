import argparse
import re
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

parser = argparse.ArgumentParser(description='select sequences with the same pattern of gaps as in the first sequence')
parser.add_argument("--input", dest='input', type=str, help="path to the file")
parser.add_argument("--output", dest='output', type=str, help="path to the file")
args = parser.parse_args()

seqs = SeqIO.parse(args.input, "fasta")

newsrs = []
for s in seqs:
	newseq = re.split("[UME]", str(s.seq))[0]
	newsr = SeqRecord(Seq(newseq, generic_dna), id=s.id, description="")
	newsrs.append(newsr)
		

SeqIO.write(newsrs, args.output, "fasta")