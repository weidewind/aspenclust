#!/usr/bin/env Rscript

library(optparse)


option_list = list(
  make_option(c("--folder"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--file"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--output"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--zero"), type="character", default=0, 
              help="sub zero p-values by this value", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
zero = as.numeric(opt$zero)

tabl = read.table(paste(c(opt$folder, "/", opt$file), collapse = ""), header = T, sep = " ", stringsAsFactors = F)

subzero = function(vec, newzero){
	res = sapply(vec, function(v){
		if(v == 0){
			newzero
		}else {v}
	})
	return(res)
}

print(zero)
tabl$upval = subzero(tabl$upval, zero)
tabl$upval_adjusted = p.adjust(tabl$upval, method = "BH")
tabl$lpval = subzero(tabl$lpval, zero)
tabl$lpval_adjusted = p.adjust(tabl$lpval, method = "BH")
print("with BH")
print(tabl)
write.table(tabl, opt$output, row.names = F, quote = F)