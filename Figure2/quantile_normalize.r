library(preprocessCore)
Args <- commandArgs()
celltype = Args[6]
QN_SDOC <- function(file_input){
	data = read.table(file = paste0("./1-tmp/",file_input), sep = "\t")

	d = as.matrix(as.numeric(data[,1]))
	normalize.quantiles.use.target(d, target = rnorm(10000), copy=FALSE, subset=NULL)
	d <- as.character(d)
	data[,1] = d

	data <- as.vector(data)
	dir.create("./1-QNed_SDOC")
	write.table(file = paste0("./1-QNed_SDOC/",file_input), data, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
	return(0)
}


QN_SDOC(paste0(celltype,"_10kb_log2n_SDOCs"))
QN_SDOC(paste0(celltype,"_10kb_SDOCs"))
