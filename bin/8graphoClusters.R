rm(list = ls())
setwd("/home/clovis/Dropbox/Chumbo/")
load("./transcriptograms/allTranscriptogramers125")
load("./Data/tfs.RData")
#load("./Data/colors.RData")
load("./Data/colorsSB.RData")

figuras="figurasSB"

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

#load(file='clusters12.RData')
load(file='clusters12SB.RData')

#intersecção dos clusteres
c1<-clusteres$T1
c2<-clusteres$T2

#cria dtf com intersecção
cltRef<-data.frame(c1=c1,c2=c2)

#load(file = "./Data/cl12.RData")
load(file = "./Data/cl12SB.RData")

#extrai genes para grupos 3 4 5 com base no 3 ref3
#inicioLoop ----
ref2=2

clust<-list()
interval=1
cluster=1
for(interval in 1:2){
  for(cluster in 1:length(color[[interval]])){
    resetd(rdp1)
    
    rdp1<-clusterVisualization(transc[[interval]],onlyGenesInDE = F, 
                               clusters = c(cluster),connected = T)
    Sys.sleep(5)
    
    #RedeR::calld(rdp1)
    
    clust[[cluster]]<-getGraph(rdp1)
    
    #resetd(rdp1)
    
    edges<-data.frame(unique(get.edgelist(clust[[cluster]])),stringsAsFactors = F)
    
    colnames(edges)<-c("source","target")
    edges$interaction<-rep("interacts",nrow(edges))
    
    nodes<-data.frame(id=V(clust[[cluster]])$name, nada=rep("N",length(V(clust[[cluster]])$name)))
    # optional
    
    RCy3::deleteAllNetworks()
    RCy3::setNodeShapeDefault("ELLIPSE")
    #substr remove os parametros de alfa da cor
    cor<-substr(x = color[[interval]][cluster],start = 1,stop = 7)
    RCy3::setNodeBorderColorDefault("#AAAAAA")
    RCy3::setNodeColorDefault(substr(x = color[[interval]][cluster],start = 1,stop = 7))
    createNetworkFromDataFrames(nodes,edges)
    #RCy3::setNodeLabelMapping(table.column = "names")
    #Sem opção seta o prefered layout
    RCy3::layoutNetwork()
    RCy3::exportImage(filename = paste0("/home/clovis/Dropbox/Chumbo/",figuras,"redesFig1/itv",
                                        interval,"cl",cluster,".svg"),type = "SVG")
    # gD.cyt <- igraph::as_graphnel(clust[[i]])
    # 
    # gDCW <- RCy3::createNetworkFromGraph(title = "Les Miserables", graph = gD.cyt)
    # 
    
  }
}
#save(clust,file = "ClustersTransc3.RData")

