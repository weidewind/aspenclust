#!/usr/bin/env Rscript

library(optparse)


option_list = list(
  make_option(c("--folder_usher"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--folder_iqtree"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--cfile"), type="character", default=NULL, 
              help="${cfile}_with_disttest.csv", metavar="character"),
  make_option(c("--dfile"), type="character", default=NULL, 
              help="${dfile}_with_disttest.csv", metavar="character"),
  make_option(c("--output"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--zero"), type="character", default=0, 
              help="sub zero p-values by this value", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
zero = as.numeric(opt$zero)

iqtree_c = read.table(paste(c(opt$folder_iqtree, "/", opt$cfile), collapse = ""), header = T, sep = " ", stringsAsFactors = F)
print(iqtree_c)
iqtree_c$type = "c" #rep("c", nrow(iqtree_c))
iqtree_d = read.table(paste(c(opt$folder_iqtree, "/", opt$dfile), collapse = ""), header = T, sep = " ", stringsAsFactors = F)
print(iqtree_d)
iqtree_d$type = "d" #rep("d", nrow(iqtree_d))
print(iqtree_d)
iqtree = rbind(iqtree_c, iqtree_d)
print(iqtree)

usher_c = read.table(paste(c(opt$folder_usher, "/", opt$cfile), collapse = ""), header = T, sep = " ", stringsAsFactors = F)
usher_c$type =  "c" #rep("c", nrow(usher_c))
print(usher_c)
usher_d = read.table(paste(c(opt$folder_usher, "/", opt$dfile), collapse = ""), header = T, sep = " ", stringsAsFactors = F)
usher_d$type = "d" #rep("d", nrow(usher_d))
print(usher_d)
usher = rbind(usher_c, usher_d)
print(usher)

subzero = function(vec, newzero){
	res = sapply(vec, function(v){
		if(v == 0){
			newzero
		}else {v}
	})
	return(res)
}

print(zero)
iqtree$lpval = subzero(iqtree$lpval, zero)
iqtree$upval = subzero(iqtree$upval, zero)
iqtree$lpval_adjusted = p.adjust(iqtree$lpval, method = "BH")
iqtree$upval_adjusted = p.adjust(iqtree$upval, method = "BH")
print("with BH")
print(iqtree)
usher$lpval = subzero(usher$lpval, zero)
usher$upval = subzero(usher$upval, zero)
usher$lpval_adjusted = p.adjust(usher$lpval, method = "BH")
usher$upval_adjusted = p.adjust(usher$upval, method = "BH")
merged = merge(iqtree, usher, by = c("site1", "site2"), all.x = T, all.y = T)
write.table(merged, opt$output, row.names = F, quote = F)