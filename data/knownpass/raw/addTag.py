import argparse
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

parser = argparse.ArgumentParser(description='select sequences with the same pattern of gaps as in the first sequence')
parser.add_argument("--passages", dest='passages', type=str, help="path to the file")
parser.add_argument("--input", dest='input', type=str, help="path to the file")
parser.add_argument("--output", dest='output', type=str, help="path to the file")
args = parser.parse_args()

with open(args.passages, "r") as passages:
	passDict = dict([s.strip() for s in str.split("\t")] for str in passages)

align = AlignIO.read(args.input, "fasta")

newseqs = []
for s in align:
	t = passDict[s.id]
	if t == "SIAT":
		t = "UNPASSAGED"
	tag = t*3
	news = SeqRecord(Seq(str(s.seq) + tag, generic_dna), id=s.id, description="")
	newseqs.append(news)
		

SeqIO.write(newseqs, args.output, "fasta")