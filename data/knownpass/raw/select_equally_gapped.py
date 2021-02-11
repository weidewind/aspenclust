from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import argparse

parser = argparse.ArgumentParser(description='select sequences with the same pattern of gaps as in the first sequence')
parser.add_argument("--input", dest='input', type=str, help="path to the file")
parser.add_argument("--output", dest='output', type=str, help="path to the file")
parser.add_argument("--leftovers", dest='leftovers', type=str, help="path to the file")
args = parser.parse_args()

ali = AlignIO.read(args.input, "fasta")
refseq = ali[0]
refgaps=[ind for ind, c in enumerate(refseq.seq) if c == "-"]

validSeqs=[]
invalidSeqs=[]
refgapsExt=[0,1,2]
refgapsExt.extend(refgaps)

for a in ali[1:]:
	agaps=[ind for ind, c in enumerate(a.seq) if c == "-"]
	if agaps == refgaps:
		validSeqs.append(a)
	elif agaps == refgapsExt:
		atga=a
		atga.seq="atg"+a.seq[3:]
		validSeqs.append(atga)
	else:
		invalidSeqs.append(a)

with open(args.output, "w") as out:
	align = MultipleSeqAlignment(validSeqs)
	AlignIO.write(align, out, "fasta")
with open(args.leftovers, "w") as out:
	align = MultipleSeqAlignment(invalidSeqs)
	AlignIO.write(align, out, "fasta")