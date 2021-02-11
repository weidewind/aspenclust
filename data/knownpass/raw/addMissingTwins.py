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
parser.add_argument("--iqlog", dest='iqlog', type=str, help="path to the file")
parser.add_argument("--output", dest='output', type=str, help="path to the file")
args = parser.parse_args()

align = AlignIO.read(args.fasta, "fasta")
alignIndex = dict((s.id, ind) for ind, s in enumerate(align))
validSeqs = []

with open(args.iqlog, "r") as missing:
    for s in [x.strip() for x in missing]:
        if s[-16:] == "added at the end":
            missid, presentid = s.split()[1:5:3]
            presentid = presentid[:-1]
            print(presentid)
            print(alignIndex[presentid])
            seq = align[alignIndex[presentid]].seq
            misseq = SeqRecord(seq, id=missid, description="")
            align.append(misseq)
            print("missing " + missid + " added")

with open(args.output, "w") as output:
    newalign = MultipleSeqAlignment(align)
    AlignIO.write(newalign, output, "fasta")