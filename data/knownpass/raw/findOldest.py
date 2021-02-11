import re
import argparse
from datetime import datetime
from Bio import AlignIO

parser = argparse.ArgumentParser()
parser.add_argument("--datefile", dest='datefile', type=str, help="path to the all.fa file")
parser.add_argument("--fasta", dest='fasta', type=str, help="path to the file with passages")
args = parser.parse_args()


def toDate(datestr):
    if not re.match("^[0-9-]+$", datestr):
        halfdate = datestr.split("_")[0]
        if "-" in halfdate:
            datestr = halfdate + "-28"
        else:
            datestr = halfdate + "-12-31"
    try:
        mydate = datetime.strptime(datestr, "%Y-%m-%d")
    except ValueError as v:
        print(v)
        print(datestr)
        mydate = datetime.strptime("3000-01-01", "%Y-%m-%d")
    return (mydate)


dateDict = {}

with open(args.datefile) as datef:
    for line in datef:
        line = line.strip()
        if line.startswith(">"):
            year = re.split('[-_]', line.split("=")[5])[0]
            idENTs = list(filter(lambda x: 'ENTIFIER' in x, line.split("=")))
            if len(idENTs) > 0:
                id = "STRAIN" + str(idENTs[0][:-8]) + "_" + year
                fulldate = line.split("=")[5]
                dateDict[id] = toDate(fulldate)

ali = AlignIO.read(args.fasta, "fasta")
print(sorted(ali, key=lambda a: dateDict[a.id])[0].id)

#print(sorted(dateDict.items(), key=lambda i: i[1])[:10])
