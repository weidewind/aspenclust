import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", dest='input', type=str, help="path to the all.fa file")
parser.add_argument("--output", dest='output', type=str, help="path to the file with passages")
args = parser.parse_args()

with open(args.input) as inf:
	seqs = {}
	for line in inf:
		if not line.startswith("#"):
			id, site, state = line.split()[:3]
			if not id == "Node":
				if id not in seqs:
					seqs[id] = ""
				seqs[id] = seqs[id] + state

with open(args.output, "w") as output:
	for item in seqs.items():
		output.write(">" + item[0] + "\n")
		output.write(item[1] + "\n")
