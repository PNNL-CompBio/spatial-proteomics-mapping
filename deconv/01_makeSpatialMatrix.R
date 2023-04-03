##loads spatial data

setwd("../")
source("spleenDataFormatting.R")
setwd("deconv")

##write out spatial data, with exponential call
mat<-exprs(spat.prot)
expmat<-2^mat
write.table(as.matrix(mat),file='spatialProtMat.tsv',sep='\t',row.names=T,col.names=T,quote=F)

##write out sorted data
mat<-exprs(global.sorted)
expmat<-2^mat
write.table(as.matrix(mat),file='globalProtMat.tsv',sep='\t',row.names=T,col.names=T,quote=F)