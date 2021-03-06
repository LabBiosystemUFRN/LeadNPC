rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadNPC/"
setwd(baseDir)
source(file = "bin/00base.R")

load("./Data/allTranscriptogramers80")
load("./Data/tfs.RData")
load("./Data/colors.RData")

figures="tmpGraf"

library(igraph)
library(transcriptogramer)
library(ggplot2)
library(purrr)
library(scales)
library("RCy3")



#intersecção dos clusteres
load(file = "./Data/clusteres.RData")
c1<-clusteres$T1
c2<-clusteres$T2

load(file='./Data/clusters12.RData')

#intersecção dos clusteres
c1<-clusteres$T1
c2<-clusteres$T2

#cria dtf com intersecção
cltRef<-data.frame(c1=c1,c2=c2)

load(file = "./Data/cl12.RData")
#load(file = "./Data/cl12SB.RData")

#extrai genes para grupos 3 4 5 com base no 3 ref3
#inicioLoop ----
ref2=2

clust<-list()
interval=1
cluster=3
assoc<-transc[[1]]@association

for(interval in 1:2){
  limite<-length(color[[interval]])-1 # com boundary
#  limite<-length(color[[interval]])
  for(cluster in 1:limite){
    options(stringsAsFactors = FALSE)
    #cria dataframe com nodes e arestas    
    if(interval==1){
      nodes<-subset(clusters12,Clust1==cluster,select = "Protein") 
    }else{
      nodes<-subset(clusters12,Clust2==cluster,select = "Protein") 
    }
    nodes$Protein<-as.character(nodes$Protein)
    colnames(nodes)<-c("id")
    nodes$name<-rep(".",nrow(nodes))
    
    edges<-subset(assoc, p1%in%nodes$id & p2%in%nodes$id)
    edges<- data.frame(t(apply(edges,1,sort)))
    edges<-unique(edges)
    
    colnames(edges)<-c("source","target")
    edges$interaction<-rep("interacts",nrow(edges))
    
    
    RCy3::deleteAllNetworks()
    RCy3::setNodeShapeDefault("ELLIPSE")
    #substr remove os parametros de alfa da cor
    cor<-substr(x = color[[interval]][cluster],start = 1,stop = 7)
    createNetworkFromDataFrames(nodes,edges)
    RCy3::setNodeBorderColorDefault("#888888")
    RCy3::setNodeBorderWidthDefault(8)
    RCy3::setNodeLabelColorDefault(cor)
    #Sem opção seta o prefered layout
    RCy3::layoutNetwork()
    RCy3::setNodeColorDefault(cor,style.name = "Chumbo")
    RCy3::setEdgeColorDefault(cor,style.name = "Chumbo")
    RCy3::setEdgeLineWidthDefault(10,style.name = "Chumbo")
    RCy3::setVisualStyle("Chumbo")
    RCy3::exportImage(filename = paste0("./",figures,"/t",
                                        interval,"c",cluster,".svg"),type = "SVG")

  }
}
#save(clust,file = "ClustersTransc3.RData")

