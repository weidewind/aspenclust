import argparse
from Bio import AlignIO

parser = argparse.ArgumentParser(description='select sequences with the same pattern of gaps as in the first sequence')
parser.add_argument("--fasta", dest='fasta', type=str, help="path to the file")
parser.add_argument("--passages", dest='passages', type=str, help="path to the file")
parser.add_argument("--output", dest='output', type=str, help="path to the file")
args = parser.parse_args()

with open(args.passages, "r") as passages:
	passDict = dict([s.strip() for s in str.split("\t")] for str in passages)

align = AlignIO.read(args.fasta, "fasta")
passCode = {"SIAT": 0, "UNPASSAGED": 0, "MDCK": 1, "EGG": 2, "MONKEYKIDNEY": 3}

with open(args.output, "w") as output:
	for a in align:
		output.write(">" + a.id + "\n")
		output.write(str(passCode[passDict[a.id]]) + "\n")
