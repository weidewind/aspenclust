library(ape)
library(plyr)
library(stringr)

isdate<- function(date){
  if (str_detect(date, "^[0-9]{4}-[0-9]{2}-[0-9]{2}$")){
    TRUE
  }
  else{
    FALSE
  }
}


treefile = "data/noeggs/n1pandemic.l.r.newick"
datefile = "../evopoisson/data/noeggs_ksu/raw_and_processed/fulldaten1pdm_passages"

min.date = 14300 #circa 2009
nsteps = 20
output = paste(c("data/noeggs/","n1pandemic_noeggs_dates_nsteps", nsteps),collapse = "")

tree<-read.tree(file = treefile, text = NULL, tree.names = NULL, skip = 0,
         comment.char = "", keep.multi = FALSE)
dates<-read.csv2(file=datefile, sep = "\t", header = FALSE,stringsAsFactors = FALSE)


strain_date<-ldply(dates[[1]], function(g){
  splr <-strsplit(g,"_")[[1]]
  row<-data.frame(strain = as.character(splr[1]), date = as.character(splr[2]), stringsAsFactors = FALSE)
})

node.dates <-sapply(tree$tip.label, function(e){
  strain = strsplit(e, "_")[[1]][1]
  date <- strain_date[strain_date$strain == strain,2][[1]]
  if (!isdate(date)){
    regex = paste(c("^", strain, "$"), collapse = "")
    rownum = grep(regex,strain_date$strain)
    while (!isdate(date)){
      rownum = rownum-1
      date = strain_date[rownum,2][[1]]
    }
  }
  as.Date(date) # stored internally as the number of days since January 1, 1970
})

mu = estimate.mu(tree, node.dates, p.tol = 0.05)
is.binary = FALSE #is.binary = is.binary.phylo(tree)
estimation<-estimate.dates(tree, node.dates, mu = mu,
               min.date = min.date, show.steps = 1,
               opt.tol = 1e-8, lik.tol = 0,
               nsteps = nsteps,is.binary = is.binary)

tree$edge.length <- estimation[tree$edge[, 2]]
write.tree(tree, output)




