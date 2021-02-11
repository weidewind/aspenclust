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
description='select sequences with the same pattern of gaps as in the first sequence'
    )
parser.add_argument("--tree", dest='tree', type=str, help="path to the file")
parser.add_argument("--phylip", dest='phylip', type=str, help="path to the file")
parser.add_argument("--fasta", dest='fasta', type=str, help="path to the file")
parser.add_argument("--output", dest='output', type=str, help="path to the file")
args = parser.parse_args()

tree = Tree(args.tree, format = 1)
if (args.phylip):
        align = AlignIO.read(args.phylip, "phylip-relaxed")
else:
        align = AlignIO.read(args.fasta, "fasta")
alignIndex = dict((s.id, ind) for ind, s in enumerate(align))

ambiDict = {}
for seq in align:
    inds = []
    for ind, c in enumerate(seq.seq):
        if c not in "atgcATGC":
            inds.append(ind)
    if len(inds) > 0:
        ambiDict[seq.id] = inds

disambiDict = collections.defaultdict(dict)

# def goneTooFar(tnode):
#     print (mindist)
#     if (tnode.get_distance(node) > mindist):
#         return(True)
#     else:
#         return(False)

def myTraverse(node, curdist, mindist, newchar, ind): #curdepth - current distance to the leaf with ambichar; mindist - best distance known as far
   # print (" ".join(["curdist", str(curdist), "mindist", str(mindist),"newchar", newchar, "ind", str(ind)]))
    if node.is_leaf() and curdist < mindist:
        tchar = align[alignIndex[node.name]][ind]
        if tchar in "atgcATGC":
            mindist = curdist
            newchar = tchar
    elif curdist < mindist:
        for ch in node.children:
            chdist = curdist + ch.dist
            if chdist < mindist:
                if ch.is_leaf():
                    tchar = align[alignIndex[ch.name]][ind]
                    if tchar in "atgcATGC":
                        mindist = chdist
                        newchar = tchar
            #            return(mindist, newchar) #nope, out!
            #        else:
            #            return([])    #nope, out!
                else:
                    out = myTraverse(ch, chdist, mindist, newchar, ind)
                    if out[0] < mindist and out[1] in "atgcATGC":
                        mindist = out[0]
                        newchar = out[1]
           # print (" ".join(["node in question: ", ch.name, "chdist", str(chdist), "mindist", str(mindist),"newchar", newchar, "ind", str(ind)]))
    return(mindist, newchar)

  


for seqid in ambiDict.keys():
    print("seqid " + seqid)
    node = tree.get_leaves_by_name(name=seqid)[0]
    for ind in ambiDict[seqid]:
        tnode = node
        mindist = 1000.0
        newchar = "X"
        fixed = False
        while(not tnode.is_root() and (not fixed or node.get_distance(tnode) < mindist)):
            siss = tnode.get_sisters()
            for s in siss:
              #  print ("sis "+s.name)
                curdist = node.get_distance(s)
                sismindist, sisnewchar = myTraverse(s, curdist, mindist, newchar, ind)
                if sismindist < mindist:
                    mindist = sismindist
                    newchar = sisnewchar
            if newchar in "atgcATGC":
                fixed = True
            tnode = tnode.up
      #  print("seq "+seqid+" ind "+str(ind)+" node "+node.name+" mindist "+ str(mindist)+" newchar "+newchar)    
        disambiDict[seqid][ind] = newchar

newseqs = []
for a in align:
    if a.id in disambiDict.keys():
        newseq = str(a.seq)
        for ind in ambiDict[a.id]:
            newseq = newseq[:ind] + disambiDict[a.id][ind] + newseq[ind+1:]
        newa = SeqRecord(Seq(newseq, generic_dna), id=a.id, description="")
        newseqs.append(newa)
    else:
        newseqs.append(a)

with open(args.output, "w") as output:
    newalign = MultipleSeqAlignment(newseqs)
    AlignIO.write(newalign, output, "fasta")
