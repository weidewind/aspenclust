import argparse
import re
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
import collections

# node = None
# mindist = 1000

parser = argparse.ArgumentParser(
description='select sequences with the same pattern of gaps as in the first sequence')

parser.add_argument("--fasta", dest='fasta', type=str, help="path to the file")
parser.add_argument("--output", dest='output', type=str, help="path to the file")
args = parser.parse_args()

align = AlignIO.read(args.fasta, "fasta")
validSeqs = [s for s in align if len(s.translate(to_stop=True)) == len(s)/3]

with open(args.output, "w") as output:
	newalign = MultipleSeqAlignment(validSeqs)
	AlignIO.write(newalign, output, "fasta")