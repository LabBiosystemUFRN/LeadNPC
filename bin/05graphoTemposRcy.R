rm(list = ls())

figuras="figurasSB"

setwd("/home/clovis/Dropbox/Chumbo/")
load("./transcriptograms/allTranscriptogramers80")
load("./Data/tfs.RData")
#load("./Data/colors.RData")
load("./Data/colorsSB.RData")

library(igraph)
library(RedeR)
library(transcriptogramer)
library(ggplot2)
library(purrr)
library(scales)
library("RCy3")


myTheme<-list()
zoom=30
myTheme$nestShape='ROUNDED_RECTANGLE'
myTheme$nestColor='#ffffff'
myTheme$nestLineWidth=4*(100/zoom)
myTheme$nestFontSize=24*(100/zoom)
myTheme$nestFontX=5
myTheme$nestFontY=10.8
myTheme$isAssign=TRUE


#intersecção dos clusteres
#load(file = "./Data/clusteres.RData")
load(file = "./Data/clusteresSB.RData")
c1<-clusteres$T1
c2<-clusteres$T2

#load(file='./Data/clusters12.RData')
load(file='./Data/clusters12SB.RData')

#intersecção dos clusteres
c1<-clusteres$T1
c2<-clusteres$T2

#cria dtf com intersecção
cltRef<-data.frame(c1=c1,c2=c2)


load(file = "./Data/cl12.RData")
load(file = "./Data/cl12SB.RData")

#extrai genes para grupos 3 4 5 com base no 3 ref3
#inicioLoop ----
ref2=2
clust<-list()
interval=1
cluster=5
assoc<-transc[[1]]@association
for(interval in 1:2){
  for(cluster in 1:length(color[[interval]])){
    interval=2
    cluster=4
    options(stringsAsFactors = FALSE)
    
    # rdp1 <- RedPort()
    # resetd(rdp1)
    # 
    # rdp1<-clusterVisualization(transc[[interval]],onlyGenesInDE = F,
    #                            clusters = c(cluster),connected = T)
    # 
    # clustF<-RedeR::getGraph(rdp1)
    # 
    # n1<-V(clustF)$name
    # e1<-data.frame(unique(get.edgelist(clustF)),stringsAsFactors = F)
    
    nodesF<-subset(clusters12,Clust1==5,select = "Protein") 
    nodesF$Protein<-as.character(nodesF$Protein)
    colnames(nodesF)<-c("id")
    
    edgesF<-subset(assoc, p1%in%nodesF$id & p2%in%nodesF$id)
    edgesF<- data.frame(t(apply(edgesF,1,sort)))
    edgesF<-unique(edgesF)
    
    colnames(edgesF)<-c("source","target")
    edgesF$interaction<-rep("interacts",nrow(edgesF))
    
    
    #Pai
    nodesP<-subset(clusters12,Clust2==4,select = c("Protein") )
    nodesP$Protein<-as.character(nodesP$Protein)
    colnames(nodesP)<-c("id")
    
    edgesP<-subset(assoc, p1%in%nodesP$id & p2%in%nodesP$id)
    edgesP<- data.frame(t(apply(edgesP,1,sort)))
    edgesP<-unique(edgesP)
    
    colnames(edgesP)<-c("source","target")
    edgesP$interaction<-rep("interacts",nrow(edgesP))
  

    RCy3::deleteAllNetworks()
    RCy3::setNodeShapeDefault("ELLIPSE")
    #substr remove os parametros de alfa da cor
    cor<-substr(x = color[[interval]][cluster],start = 1,stop = 7)
    RCy3::setNodeBorderColorDefault("#AAAAAA")
    RCy3::setNodeColorDefault(substr(x = color[[interval]][cluster],start = 1,stop = 7))
    createNetworkFromDataFrames(nodesP,edgesP)
    createGroup(group.name = "teste",nodes = c(nodesF$id),nodes.by.col = "id")
    
    createSubnetwork(nodes = c(nodesF$id),nodes.by.col = "id")
    
    
    
    #RCy3::setNodeLabelMapping(table.column = "names")
    #Sem opção seta o prefered layout
    RCy3::layoutNetwork()
    RCy3::exportImage(filename = paste0("/home/clovis/Dropbox/Chumbo/",figuras,"/redesFig1/itv",
                                        interval,"cl",cluster,".svg"),type = "SVG")
    # gD.cyt <- igraph::as_graphnel(clust[[i]])
    # 
    # gDCW <- RCy3::createNetworkFromGraph(title = "Les Miserables", graph = gD.cyt)
    # 
    
  }
}
#save(clust,file = "ClustersTransc3.RData")

