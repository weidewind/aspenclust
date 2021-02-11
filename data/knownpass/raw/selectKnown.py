import argparse
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

parser = argparse.ArgumentParser(description='select sequences with the same pattern of gaps as in the first sequence')
parser.add_argument("--fasta", dest='fasta', type=str, help="path to the file")
parser.add_argument("--passages", dest='passages', type=str, help="path to the file")
parser.add_argument("--output", dest='output', type=str, help="path to the file")
args = parser.parse_args()

with open(args.passages, "r") as passages:
	passDict = dict([s.strip() for s in str.split("\t")] for str in passages)

align = AlignIO.read(args.fasta, "fasta")

validSeqs = [s for s in align if s.id in passDict.keys() and passDict[s.id] in ["SIAT", "UNPASSAGED", "EGG", "MDCK", "MONKEYKIDNEY"]]

with open(args.output, "w") as output:
	newalign = MultipleSeqAlignment(validSeqs)
	AlignIO.write(newalign, output, "fasta")	