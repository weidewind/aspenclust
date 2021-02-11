from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
align = AlignIO.read("n1.uniq", "fasta")
print(align["STRAIN152795_2007",:])