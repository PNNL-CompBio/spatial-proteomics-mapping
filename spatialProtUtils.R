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


#' expToLongForm: helper function that moves expression
#' matrix to long form
expToLongForm<-function(sce,rowname='prot'){
  
  exprs(spat.phos)%>%
    as.matrix()%>%
    as.data.frame()%>%
    tibble::rownames_to_column(rowname)%>%
    tidyr::pivot_longer(c(2:(1+ncol(exprs(spat.phos)))),names_to='Voxel',values_to='LogRatio')
}

#' calcCorrelationWithScore: calculates the correlation of each
#' element with a numeric vector of the score, such as distance or 
#' an immune score, puts value in RowData
#' @param sce
#' @param scoreName: name of score
#' @param protVal: name of feature, e.g. prot or substrate
#' @param method: spearman or pearson
#' @return SingleCellExperiment objecti
calcCorrelationWithScore<-function(sce, 
                                   scoreName,
                                   protVal='prot',
                                   method='spearman'){
    
    #join value with expression
    score.dat<-colData(sce)%>%
      as.data.frame()%>%
      dplyr::select(all_of(scoreName))%>%
      tibble::rownames_to_column('Voxel')%>%
      left_join(expToLongForm(sce,protVal))
    
    #calculate correlation
    prot.cor<-score.dat%>%
      dplyr::rename(pname=protVal,score=scoreName)%>%
      group_by(pname)%>%
      summarize(corVal=cor(score,LogRatio,method=method))
   
    newDataFrame<-rowData(sce)
    newDataFrame[[paste(scoreName,'correlation')]]<-unlist(tibble::column_to_rownames(prot.cor,'pname'))
    rowData(sce)<-newDataFrame
    
    sce
}


#' plotFeatureGrid: plots a numeric value of a single feature or set of features
#' @param sce
#' @param features 
#' 
plotFeatureGrid<-function(sce,feats,featname){
  require(ggplot2)
  require(SingleCellExperiment)
  
  p<-calcSigScore(sce,feats,featname)%>%
    plotSigGrid(featname)

    return(p)
}


#'plotSigGrid: plots a numeric value using the Xcoord and Ycoord columns
#'takes a numeric score from the colDAta
#'@param sce SingleCellExperiment
#'
plotSigGrid<-function(sce,sigName){
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
 # res <- data.frame(featureID=rownames(res), res, stringsAsFactors = F)
  colnames(res)<-paste(paste(column,'limma'),colnames(res))
  res<-res%>%
    tibble::rownames_to_column('X')
  rd<-rowData(sce)%>%
    as.data.frame()%>%
    full_join(res)
  rowData(sce)<-rd
  return(sce)
}



#' buildnetwork
#' Builds network using SingleCellExpression object and values
#' that were applied to the row data
#' @import PCSF
buildNetwork<-function(sce,featName,addPhos=FALSE){
  library(PCSF)
  strdb<-data('STRING')
  
  ##first get diffex proteins
  
}
