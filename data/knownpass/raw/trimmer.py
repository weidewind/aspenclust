from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
h1 = AlignIO.read("h1.max2ambichar.strains.reordered.aligned", "fasta")
with open("h1.max2ambichar.strains.reordered.aligned.trimmed", "w") as f:
	validList=[s for s in h1 if s.id not in ["STRAIN29082_2000","STRAIN29215_2000","STRAIN29253_2000","STRAIN713222_2015"]] # insertions and one from 2015
	validAlign = MultipleSeqAlignment(validList)
	AlignIO.write(validAlign[:, 126:1913], f, "fasta") # as is in little_ksu

n1 = AlignIO.read("n1.max2ambichar.strains.reordered.aligned", "fasta")	
with open("n1.max2ambichar.strains.reordered.aligned.trimmed", "w") as f:
	validList=[s for s in n1 if s.id not in ["STRAIN29087_2000", "STRAIN29220_2000", "STRAIN29258_2000"]]
	validAlign = MultipleSeqAlignment(validList)
	AlignIO.write(validAlign[:, 37:1518], f, "fasta") # as is in little_ksu 
	
h3 = AlignIO.read("h3.max2ambichar.strains.reordered.aligned", "fasta")	
with open("h3.max2ambichar.strains.reordered.aligned.trimmed", "w") as f:
	validList=[s for s in h3 if s.id not in ["STRAIN324462_2011","STRAIN1122601_2012","STRAIN1041294_2016","STRAIN1230945_2014"]] # insertions and one from minus strand (h3.max2ambichar.trimmed.realigned.uniqueseq.phy produced by iqtree )
	validAlign = MultipleSeqAlignment(validList)
	AlignIO.write(validAlign[:, 152:2167], f, "fasta") # as is in little_ksu 
	
n2 = AlignIO.read("n2.max2ambichar.strains.reordered.aligned", "fasta")	
with open("n2.max2ambichar.strains.reordered.aligned.trimmed", "w") as f:
	validList=[s for s in n2 if s.id not in ["STRAIN1288169_2018"]]
	validAlign = MultipleSeqAlignment(validList)
	AlignIO.write(validAlign[:, 245:2007], f, "fasta") # as is in little_ksu 

h1pdm = AlignIO.read("h1pdm.max2ambichar.strains.reordered.aligned", "fasta")	
with open("h1pdm.max2ambichar.strains.reordered.aligned.trimmed", "w") as f:
	validList=[s for s in h1pdm if s.id not in ["STRAIN301811_2010", "STRAIN1279385_2017"]]
	validAlign = MultipleSeqAlignment(validList)
	AlignIO.write(validAlign[:, 1202:3506], f, "fasta") # as is in little_ksu 

n1pdm = AlignIO.read("n1pdm.max2ambichar.strains.reordered.aligned", "fasta")	
with open("n1pdm.max2ambichar.strains.reordered.aligned.trimmed", "w") as f:	
	validList=[s for s in n1pdm if s.id not in ["STRAIN835122_2016", "STRAIN835118_2016", "STRAIN498593_2013", "STRAIN498590_2013", "STRAIN1377834_2019"]]
	validAlign = MultipleSeqAlignment(validList)
	AlignIO.write(validAlign[:, 875:2296], f, "fasta") # as is in little_ksu h3 	

#	n1 37:1518
#	n2 245:2007
#	h1pdm 1202:3506
#	n1pdm 875:2296