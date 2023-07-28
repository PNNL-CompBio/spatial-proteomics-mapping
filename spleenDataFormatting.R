##this file formats the spleen data so it can be ingested into the format we need


library(readxl)
library(dplyr)
library(tidyr)
library(SingleCellExperiment)
library(org.Hs.eg.db)
source("synapseUtil.R")
syn<-loadSynapse()


swp.map<-as.data.frame(org.Hs.egUNIPROT)%>%
  full_join(as.data.frame(org.Hs.egSYMBOL))

global.spleen<-read.table('../../Projects/HubMap/spleen/Vocanol_tip_Global.txt',sep='\t',header=T,comment.char = '')%>%
  dplyr::select(Entry.Name,Gene,starts_with('sample'))%>%
  #subset(Gene!='')%>%##TODO FIX THIS 
  tidyr::pivot_longer(cols=starts_with('sample'),names_to='Sample',values_to='LogRatio')
  
phospho.spleen<-readxl::read_xlsx('../../Projects/HubMap/spleen/Volcaono_tip_phospho_siteAnnotated.xlsx')
  #read.table('../../Projects/HubMap/spleen/Vocanol_tip_phospho.txt',sep='\t',header=T,comment.char = '')

spleen.meta<-data.frame(Sample=c('sample.09','sample.10','sample.11','sample.12',
                                 'sample.13','sample.14','sample.15','sample.16'),
                        pulpAnnotation=c('red','red','red','red','white','white','white','white'))%>%
  mutate(source='Sorted')

md.file<-'../../Projects/HubMap/spleen/voxelMetadata.xlsx'
new.md.file<-'../../Projects/HubMap/spleen/2d manuscript/Voxel_metadata_histology.xlsx'
prot.file<-'../../Projects/HubMap/spleen/AfterBatchCorrection_global.txt'
phos.file<-'../../Projects/HubMap/spleen/Phospho_70percValidvalues.txt'


loadMetadata<-function(){
  tab<-readxl::read_xlsx(new.md.file)
  return(tab) 
}


loadCrossTab<-function(ctfile,isPhospho=FALSE){
  if(isPhospho)
    sk=1
  else
    sk=0
  cfile<-read.table(ctfile,sep='\t',header=T,skip=sk)%>%
    tidyr::pivot_longer(cols=starts_with('TMT'),names_to='Sample',values_to='logRatio')%>%
    separate(Sample,into=c('TMT','Channel Index'))%>%
    mutate(`Channel Index`=as.numeric(`Channel Index`))%>%
    mutate(`TMT set`=as.numeric(gsub('TMT','',TMT)))
  
  return(cfile)
}

meta<-loadMetadata()

prot<-loadCrossTab(prot.file)%>%
  full_join(meta,by=c('TMT set','Channel Index'))%>%
  dplyr::rename(`Protein`='T..Index')

phos<-loadCrossTab(phos.file,TRUE)%>%
  left_join(meta,by=c('TMT set','Channel Index'))%>%
  dplyr::rename(`Phosphosite`='T..Index')


library(Matrix)
library(SingleCellExperiment)

counts<-prot%>%
  dplyr::select(`Voxel Number`,Protein,logRatio)%>%
  tidyr::pivot_wider(values_from=logRatio,names_from='Voxel Number',values_fn=list(logRatio=mean))%>%
  subset(!is.na(Protein))%>%
  tibble::column_to_rownames('Protein')%>%
  as.matrix()

colData <- prot%>%
  dplyr::select(`Voxel Number`,Xcoord,Ycoord,pulpAnnotation,`TMT set`,Histology)%>%
  mutate(col=as.numeric(Xcoord),row=as.numeric(Ycoord))%>%
  distinct()%>%
  subset(!is.na(`Voxel Number`))%>%
  tibble::column_to_rownames('Voxel Number')%>%
  mutate(source='2D')

rowData <- data.frame(Protein=rownames(counts))
rownames(rowData)<-rownames(counts)

#cmat <- as.matrix(countdat[,rownames(i1_spots)])
#print(dim(cmat))
nas <- which(apply(counts,1,function(x) all(is.na(x))))
if(length(nas)>0)
  counts <- counts[-nas,]
#print(dim(cmat))
rowD <- data.frame(Protein=rowData[rownames(counts),])
colD <- colData[colnames(counts),]

#create the ojbect for the protein data
spat.prot<- SingleCellExperiment(assays=list(logcounts=as(counts, "dgCMatrix")),
                                rowData=rowD,
                                colData=colD)


pcounts<-phos%>%
  dplyr::select(`Voxel Number`,Phosphosite,logRatio)%>%
  tidyr::separate(Phosphosite,into=c('uniprot_id','mod'))%>%
  left_join(swp.map)%>%
  tidyr::unite(col='Phosphosite',symbol,mod,sep='_')%>%
  subset(!is.na(Phosphosite))%>%
  dplyr::select(-c(uniprot_id,gene_id))%>%
  tidyr::pivot_wider(values_from=logRatio,names_from='Voxel Number',values_fn=list(logRatio=mean))%>%
  tibble::column_to_rownames('Phosphosite')%>%
  as.matrix()

rowData <- data.frame(protein=rownames(pcounts))
rownames(rowData)<-rownames(pcounts)
prowD <- rowData[rownames(pcounts),]
pcolD <- colData[colnames(pcounts),]


#now create the same object for the phospho data
spat.phos<- SingleCellExperiment(assays=list(logcounts=as(pcounts, "dgCMatrix")),
                               rowData=prowD,
                               colData=pcolD)

pdat<-rowData(spat.phos)%>%
  as.data.frame()%>%
  tidyr::separate(X,into=c('Protein','Psite'),sep='_',remove=FALSE)%>%
  tidyr::unite(Protein,Psite,col='Phosphosite',sep='-')

rowData(spat.phos)<-pdat

pcounts<-global.spleen%>%
  dplyr::select(`Sample`,Gene,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from='Sample',values_fn=list(logRatio=mean))%>%
  subset(!is.na(Gene))%>%
  tibble::column_to_rownames('Gene')%>%
  as.matrix()

##nwo create the same object for global
global.sorted<-SingleCellExperiment(assays=list(logcounts=as(pcounts,'dgCMatrix')),
                                colData=tibble::column_to_rownames(spleen.meta,'Sample'),
                                rowData=rownames(pcounts))
                                
                                
##lastly let's process the full spleen phospho data
phcounts<-phospho.spleen%>%
  dplyr::select(Index,Gene,starts_with('sample'))%>%#first get the columns
  tidyr::separate(Index,into=c('prot','mod'),sep='_')%>% #separate swissprot
  tidyr::unite(col='Site',Gene,mod,sep="_")%>%
  dplyr::select(-prot)%>%
  tidyr::pivot_longer(cols=starts_with('sample'),names_to='oldsample',values_to='logRatio')

phcounts$sample<-stringr::str_replace(phcounts$oldsample,'-','.')
  
phcounts<-phcounts%>%
  dplyr::select(-oldsample)%>%
  tidyr::pivot_wider(names_from='sample',values_from='logRatio')%>%
  tibble::column_to_rownames('Site')%>%
  as.matrix()

phospho.sorted<-SingleCellExperiment(assays=list(logcounts=as(phcounts,'dgCMatrix')),
                                  colData=tibble::column_to_rownames(spleen.meta,'Sample'),
                                  rowData=rownames(phcounts))



##lastly, let's normalize the phospho values by the matchced global values


