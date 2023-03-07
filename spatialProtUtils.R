########
#spatialUtils
#
#script containing helper files for the spatialProtExperiment object
#
#


#' calcSigScore: calculates a signature score on a particular SingleCellExperiment object
#' This is a numeric score based on the mean rank of the proteins in teh signature
#' @param sce SingleCellExperiment Object
#' @param sigProts list of gene names
#' @param sigNamme name of signature
#' @import SingleCellExperiment
calcSigScore<-function(sce, sigProts,sigName){
  
  library(SingleCellExperiment)
  allranks<-apply(exprs(sce),2,function(x) rev(rank(x)))
  
  allpercs<-apply(exprs(sce),2,function(x) percent_rank(x))
  rownames(allpercs)<-rownames(exprs(sce))
  siggenes<-allpercs[intersect(rownames(allpercs),sigProts),]
  sigScore=apply(siggenes,2,mean,na.rm=TRUE)
  
  colData(sce)[[sigName]]<-sigScore
  sce
}



#'plotGrid: plots a numeric value using the Xcoord and Ycoord columns
#'@param sce SingleCellExperiment
#'
plotGrid<-function(sce,sigName){
  require(ggplot2)
  require(SingleCellExperiment)

  vars<-colData(sce)%>%
    as.data.frame()%>%
    dplyr::rename(signature=sigName)
  
  p<-ggplot(vars,aes(x=Xcoord,y=Ycoord,fill=signature))+
    geom_raster()+scale_fill_viridis_c()+
    theme_bw()+
    ggtitle(paste0(sigName,' signature'))
  
  return(p)
}


#'spatialDiffEx: does differential expression using annotations in object
#'
## here we do differential expression again
spatialDiffEx<-function(sce,column='pulpAnnotation', vals=c('red','white')){
  library(limma)
  
  #collect samples by factor
  samp1<-which(colData(sce)[[column]]==vals[1])
  samp2<-which(colData(sce)[[column]]==vals[2])
  
  fac <- factor(rep(c(2,1),c(length(samp2), length(samp1))))
  #print(fac)
  design <- model.matrix(~fac)
  #print(design)
  fit <- lmFit(exprs(sce)[,c(samp2,samp1)], design)
  fit <- eBayes(fit)
  
  # print(topTable(fit, coef=2))
  res <- topTable(fit, coef=2, number=Inf, sort.by="P")
  res <- data.frame(featureID=rownames(res), res, stringsAsFactors = F)
  return(res)
}
