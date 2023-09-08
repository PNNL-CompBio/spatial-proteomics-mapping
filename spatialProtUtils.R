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
  if(length(sigProts)==1)
    sigScore<-siggenes
  else
    sigScore=apply(siggenes,2,mean,na.rm=TRUE)
  
  colData(sce)[[sigName]]<-sigScore
  sce
}



#' expToLongForm: helper function that moves expression
#' matrix to long form
expToLongForm<-function(sce,rowname='prot',spotAnnotes=c()){
  
  res<-exprs(sce)%>%
    as.matrix()%>%
    as.data.frame()%>%
    tibble::rownames_to_column(rowname)%>%
    tidyr::pivot_longer(c(2:(1+ncol(exprs(sce)))),names_to='Voxel',values_to='LogRatio')
  
  if(length(spotAnnotes)>0){
    cd<-colData(sce)|>
      as.data.frame()|>
      tibble::rownames_to_column('Voxel')|>
      dplyr::select(c('Voxel',spotAnnotes))
    res<-res|>
      left_join(cd)
  }
  res<-subset(res,Voxel!='NA')
  
  return(res)
  
    
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
plotFeatureGrid<-function(sce,feats,featname,...){
  require(ggplot2)
  require(SingleCellExperiment)
  
  feats<-intersect(feats,rownames(sce))
  if(length(feats)==0)
    return(NULL)
  
  p<-calcSigScore(sce,feats,featname)|>
    plotSigGrid(featname,...)

    return(p)
}

#' plotAnnotationGrid: plots a categorical variable on the grid
#' @param sce
#' @param features 
#' 
plotAnnotationGrid<-function(sce,annotation,color_list){
  require(ggplot2)
  require(SingleCellExperiment)
  
  ##get the variables
  vars<-colData(sce)%>%
    as_tibble()|>
    as.data.frame()%>%
    dplyr::rename(signature=gsub('-','.',annotation))##i cant get rid of auto naming
  
  p<-ggplot(vars,aes(x=Xcoord,y=Ycoord,fill=signature))+
    geom_raster()+scale_fill_manual(values=color_list)+
  theme_bw()+
    ggtitle(paste0(annotation,' annotation'))
  
  return(p)
  
}

#'plotSigGrid: plots a numeric value using the Xcoord and Ycoord columns
#'takes a numeric score from the colDAta
#'@param sce SingleCellExperiment
#'
plotSigGrid<-function(sce,sigName,disFeature){
  require(ggplot2)
  require(SingleCellExperiment)

  vars<-colData(sce)%>%
    as_tibble()|>
    as.data.frame()%>%
    dplyr::rename(signature=gsub('-','.',sigName))##i cant get rid of auto naming
  
  if(!missing(disFeature))
    vars<-vars|>
    dplyr::rename(feat=disFeature)
  
  p<-ggplot(vars,aes(x=Xcoord,y=Ycoord,fill=signature))+
    geom_raster()+scale_fill_gradient(low='lightgrey',high='goldenrod')
    theme_bw()+
    ggtitle(paste0(sigName,' signature'))
  
  if(!missing(disFeature))
    p<-p+geom_point(aes(shape=feat,alpha=0.5),size=1)
  return(p)
}


#'spatialDiffEx: does differential expression using annotations in object
#'
## here we do differential expression again
spatialDiffEx<-function(sce,column='pulpAnnotation', vals=c('red','white'),feat='X'){
  library(limma)
  
  #collect samples by factor
  samp1<-which(colData(sce)[[column]]==vals[1])
  if(length(vals)>1){
    samp2<-which(colData(sce)[[column]]==vals[2])
  }else{
    samp2=setdiff(1:ncol(colData(sce)),samp1)
  }
  
  fac <- factor(rep(c(2,1),c(length(samp2), length(samp1))))
  if('TMT set'%in%colnames(colData(sce))){
    tmt<-factor(colData(sce)[c(samp2,samp1),'TMT set'])
    #print(fac)
    design <- model.matrix(~fac)#+tmt)
  }
  else{
    design <- model.matrix(~fac)
  }
  #print(design)
  fit <- lmFit(exprs(sce)[,c(samp2,samp1)], design)
  fit <- eBayes(fit)
  
  # print(topTable(fit, coef=2))
  res <- topTable(fit, coef=2, number=Inf, sort.by="P")
 # res <- data.frame(featureID=rownames(res), res, stringsAsFactors = F)
  colnames(res)<-paste(paste(column,'limma'),colnames(res))
  res<-res%>%
    tibble::rownames_to_column(feat)
  
  rd<-rowData(sce)%>%
    as.data.frame()%>%
    full_join(res)
  rowData(sce)<-rd
  return(sce)
}



#' buildnetwork
#' In the future we will build a network using the spatialExpression object, but
#' for now we just take mean values of proteins nad phosphosites
#' that were applied to the row data
#' @import PCSF
buildNetwork<-function(sig.prot.vals=c(),
                       sig.phos.vals=c(),
                       all.prot.vals=c(),
                       all.phos.vals=c(),
                       beta=.5,nrand=100,featName=''){
   require(dplyr)
  
  ##emtpy weight vectors
  # phos.vals<-c()
  # prot.vals<-c()
  # 
  # if('Phosphosite'%in%names(rowData(sce))){
  #   phos.vals<-rowData(sce)[,featName]
  #   names(phos.vals)<-rowData(sce)[,'Phosphosite']
  #   phos.vals<-phos.vals[which(phos.vals>0)]
  # }
  # 
  # if('Protein'%in%names(rowData(sce))){
  #   prot.vals<-rowData(sce)[,featName]
  #   names(prot.vals)<-rowData(sce)[,'Protein']
  #   prot.vals<-prot.vals[which(prot.vals>0)]
  # }
  # 
  if(!require('PCSF')){
    remotes::install_github('sgosline/PCSF')
    require('PCSF')
  }
  data("STRING")
  allvals<-c()
  if(length(sig.phos.vals)>0){
    #read kinase substrate database stored in data folder
    KSDB <- read.csv('PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',
                                  stringsAsFactors = FALSE)
  
    kdat<-KSDB%>%group_by(GENE)%>%
      dplyr::select(SUB_GENE,SUB_MOD_RSD)%>%
      rowwise()%>%
      dplyr::mutate(subval=paste(SUB_GENE,SUB_MOD_RSD,sep='-'))
    
    allvals<-unique(kdat$subval)
    
    mval<-mean(STRING$cost)
    adf<-apply(kdat,1,function(x)
      #for each substrate interaction, add a link from the kinase gene -> substreate -> substrate gene
      data.frame(from=c(x[['GENE']],x[['subval']]),to=c(x[['subval']],x[['SUB_GENE']]),
                 cost=c(mval,mval/8)))%>%  ##arbitrary costs based on mean cost of edges around network
      do.call(rbind,.)
    
    ##new feature: if there are substrates without kinases, still add them. 
    missed<-setdiff(names(sig.phos.vals),union(adf$from,adf$to))
    print(paste('found',length(missed),'phosphosites that dont have kinases, adding to substrates'))
    missed_subs<-sapply(missed,function(x) unlist(strsplit(x,split='-'))[1])
    #print(missed_subs)
    newdf<-data.frame(from=names(missed_subs),to=missed_subs,cost=rep(mval/8,length(missed_subs)))
   # print(head(newdf))
    adf<-rbind(adf,unique(newdf))
          
  }else{
    adf<-data.frame()
  } 
  
  ##first get diffex proteins
  ppi <- construct_interactome(rbind(STRING,adf))
  
  ##now run the code
  weights=c(sig.phos.vals,sig.prot.vals)
  allvals=c(all.prot.vals,all.phos.vals)
  #print(terms)
  subnet<-NULL
  try(
    subnet <- PCSF_rand(ppi,abs(weights), n=nrand, r=0.2,w = 10, b = beta, mu = 0.0001)
  )
  
  if(is.null(subnet))
    return("")
  
  lfcs<-allvals[match(names(V(subnet)),names(allvals))]
  lfcs[is.na(lfcs)]<-0.0
  
  ##assign proteins first
  types<-rep('proteins',length(names(V(subnet))))
  
  names(types)<-names(V(subnet))
  
  ##then assign phosphosites
  types[intersect(names(V(subnet)),names(all.phos.vals))]<-'phosphosite'
  
 
  subnet<-igraph::set.vertex.attribute(subnet,'logFoldChange',value=lfcs)
  subnet<-igraph::set.vertex.attribute(subnet,'nodeType',value=types)
  subnet<-igraph::set.edge.attribute(subnet,'interactionType',value='protein-protein interaction')
  
  
  write_graph(subnet,format='gml',file=paste0(featName,'network.gml'))
  return(list(graph=subnet,fname=paste0(featName,'.networkgml')))
}



#' adjustPhosphoWithGlobal
#' For some phospho analysis, we might want to investigate the phosphosites that are changing
#' indpendently of the proteomiocs data. Therefore we consume two independent data objects and subtract one from the other
#' @param phos.obj phosphoproteomic data
#' @param prot.obj proteomic data
#' @return phos.obj
adjustPhophoWithGlobal<-function(phos.obj,prot.obj){
  
  #first match column headeers
  samps <- intersect(colnames(exprs(phos.obj)),colnames(exprs(prot.obj)))
  
  #change new expresion values to subtract protein values when available
  newExpr<-expToLongForm(phos.obj,'Phosphosite')%>%
    tidyr::separate(Phosphosite,into=c('Protein','Site'),sep='_')%>%
    dplyr::rename(oldLogRatio='LogRatio')%>%
    left_join(expToLongForm(prot.obj,'Protein'))%>%
    mutate(diff=(oldLogRatio-LogRatio))%>%
    subset(!is.na(diff))%>%
    tidyr::unite(c(Protein,Site),col='Phosphosite',sep='_')%>%
    dplyr::select(-c(oldLogRatio,LogRatio))%>%
    tidyr::pivot_wider(names_from=Voxel,values_from=diff)%>%
    tibble::column_to_rownames('Phosphosite')%>%
    as.matrix()
    
  #create new singleCellExperiment
  new.phos<-SingleCellExperiment(assays=list(logcounts=as(newExpr,'dgCMatrix')),
                                         colData=colData(phos.obj),
                                         rowData=rownames(newExpr))
  
  new.phos
  
  
}


plotResult<-function(enrich_res){
  ###the columns dpeendon the output a bit
  library(ggplot2)
  
  odds<-enrich_res|>
    arrange(desc(oddsratio))
  
  ##get the top 20
  odds<-odds[1:min(20,nrow(enrich_res)),]|>
    tibble::rownames_to_column('Pathway')|>
    mutate(logPval=(-1*log10(pvalue)))
  
  odds$Pathway<-factor(odds$Pathway,levels=rev(odds$Pathway))
  
  ##plot
  res<-ggplot(odds,aes(x=Pathway,y=oddsratio,fill=logPval))+
    geom_bar(stat='identity')+
    coord_flip()
  
  res
}


plotKsea<-function(full.res,pvalthresh=0.05){
  library(ggplot2)
  sigRes<-subset(full.res,p.value<pvalthresh)|>
    dplyr::select(Kinase.Gene,mS,p.value)|>distinct()|>
    mutate(up=(mS>0))|>
             mutate(logP=log10(p.value))
  p<-ggplot(sigRes,aes(y=mS,x=reorder(Kinase.Gene,mS),fill=logP))+geom_bar(stat='identity')+coord_flip()
  return(p)
}

##assume we have some phospho diffex values?
ksea<-function(sce,feature='Phosphosite',FC='Histology.limma.logFC',p='Histology.limma.adj.P.Val'){
  library(KSEAapp)
  ks<-readr::read_csv("https://raw.githubusercontent.com/casecpb/KSEA/master/PSP%26NetworKIN_Kinase_Substrate_Dataset_July2016.csv")
  
  phosDiff<-rowData(sce)|>as.data.frame()
  PX<-phosDiff|>
    dplyr::rename(p=p,FC=FC)|>
    dplyr::rename(feature=feature)|>
    dplyr::mutate(Peptide='')|>
    tidyr::separate(feature,into=c("Protein",'Residue.Both'),sep='_')
  
  if(!'Gene'%in%colnames(PX))
    PX$Gene<-PX$Protein
  PX<-PX|>dplyr::select(Protein,Gene,Peptide,Residue.Both,p,FC)
  
  KSEA.Complete(ks, PX, NetworKIN=TRUE, NetworKIN.cutoff=1, m.cutoff=1, p.cutoff=0.05)
  res<-readr::read_csv('KSEA Kinase Scores.csv')
  links<-readr::read_csv("Kinase-Substrate Links.csv")
  full.res<-res|>left_join(links,by='Kinase.Gene')
  
  return(full.res)
}

gsea<-function(sce,gl,feature="Protein",FC='Histology.limma.logFC'){
  library(GSVA)
  
}
