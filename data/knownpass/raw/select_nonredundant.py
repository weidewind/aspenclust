import argparse
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

parser = argparse.ArgumentParser(description='select sequences with the same pattern of gaps as in the first sequence')
parser.add_argument("--clusters", dest='clusters', type=str, help="path to the file")
parser.add_argument("--passages", dest='passages', type=str, help="path to the file")
parser.add_argument("--input", dest='input', type=str, help="path to the file")
parser.add_argument("--output", dest='output', type=str, help="path to the file")
args = parser.parse_args()

with open(args.passages, "r") as passages:
	passDict = dict([s.strip() for s in str.split("\t")] for str in passages)

uniqIDs = []
isfirst = 1
passList = []
with open(args.clusters, "r") as clusters:
	for str in clusters:
		if (len(str) == 0):
			continue
		if (str.startswith(">")):
			print(">")
			isfirst = 1
			passList.clear()
		else:
			id = str.split()[2][1:-3]
			if (isfirst == 1):
				uniqIDs.append(id)
				passList.append(passDict[id])
				isfirst = 0
			elif (passDict[id] not in passList):
				print(passDict[id]+" not in "+",".join(passList))
				uniqIDs.append(id)
				passList.append(passDict[id])

ali = AlignIO.read(args.input, "fasta")
validSeqs = []
for a in ali:
	if a.id in uniqIDs:
		validSeqs.append(a)
with open(args.output, "w") as output:
	align = MultipleSeqAlignment(validSeqs)
	AlignIO.write(align, output, "fasta")
