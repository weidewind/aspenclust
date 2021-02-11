#!/usr/bin/env Rscript
libloc = file.path("~", "R","library", fsep = .Platform$file.sep)
list.of.packages <- c("ape", "optparse")
new.packages <- setdiff(list.of.packages, installed.packages(lib.loc = libloc)[,"Package"])
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', lib = libloc)
#install.packages("ape")
#install.packages("optparse")
library(ape, lib.loc = libloc)
library(optparse, lib.loc = libloc)


option_list = list(
  make_option(c("-t", "--treefile"), type="character", default=NULL, 
              help="protein: h1, h3, n1 or n2", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="initialization method: clusterization of parameters (cluster) or radomly chosen parameters (random) [default= %default]", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

with (opt,{
treefile <-"C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/KseniaTrees/little_ksu/n1pandemic.l.r.newick"
output<- "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/KseniaTrees/little_ksu/n1pandemic_distance_matrix.csv"
print (treefile)
print (output)
tree<-read.tree(treefile)
#PatristicDistMatrix<-cophenetic(tree) #between leafs
PatristicDistMatrix<-dist.nodes(tree) # between all nodes
#hist(PatristicDistMatrix[lower.tri(PatristicDistMatrix)],breaks=seq(from=0, to=400, by=10))
dimnames(PatristicDistMatrix) = list(c(tree$tip.label, tree$node.label), c(tree$tip.label, tree$node.label))
write.csv(PatristicDistMatrix, output)
})
