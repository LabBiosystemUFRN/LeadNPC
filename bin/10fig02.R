rm(list = ls())

setwd("/home/clovis/Dropbox/Chumbo/")
load("./transcriptograms/allTranscriptogramers80")
#load(file = "./Data/associationHs700.RData")
load("./Data/colors.RData")

DE1<-transc[[1]]@DE
DE2<-transc[[2]]@DE
setwd("/home/clovis/Dropbox/Chumbo/figuras/Figura2/redes/")
nodes<-read.table(file = "ndt2c1.txt",sep="\t",header = T, stringsAsFactors = F)
#normaliza coordenadas dos nodes entre 50 e 974
nodes$x<-(nodes$x-min(nodes$x))/(max(nodes$x)-min(nodes$x))*974+50
nodes$y<-(nodes$y-min(nodes$y))/(max(nodes$y)-min(nodes$y))*974+50

nodes$Protein<-substr(nodes$alias,1,20)
nodes$net<-substr(nodes$alias,24,25)
nodes$alias<-NULL
nodes$color[nodes$net==0]<-toupper(color[[2]][1])
nodes$color[nodes$net==1]<-toupper(color[[1]][1])
nodes$color[nodes$net==2]<-toupper(color[[1]][2])
nodes$color[nodes$net==3]<-toupper(color[[1]][3])

#expressão diferencial
# nodes1<-merge(nodes,DE1[,c(1,4,7)], by="Protein")
# nodes2<-merge(nodes,DE2[,c(1,4,7)], by="Protein")
#todos os nós
nodes1<-merge(nodes,DE1[,c(1,4,7)], by="Protein")
nodes2<-merge(nodes,DE2[,c(1,4,7)], by="Protein")


#normaliza valores de expressão entre 0 e 1
minz<-min(min(nodes1$logFC),min(nodes2$logFC))
maxz<-max(max(nodes1$logFC),max(nodes2$logFC))
nodes1$logFC<-(nodes1$logFC-minz)/(maxz-minz)
nodes2$logFC<-(nodes2$logFC-minz)/(maxz-minz)


edges<-read.table(file = "edt2c1.txt",sep="\t",header = T, stringsAsFactors = F)
edges1<-edges[edges$node_a%in%nodes1$node_id &
              edges$node_b%in%nodes1$node_id,]
edges1<-merge(edges1,nodes1[,c(2,6,7,10)], by.x = "node_a", by.y = "node_id")
edges1<-merge(edges1,nodes1[,c(2,6,7,10,5)], by.x = "node_b", by.y = "node_id")
edges1$node_b<-NULL
edges1$node_a<-NULL
edges1$weight<-NULL

edges2<-edges[edges$node_a%in%nodes2$node_id &
                edges$node_b%in%nodes2$node_id,]
edges2<-merge(edges2,nodes2[,c(2,6,7,10)], by.x = "node_a", by.y = "node_id")
edges2<-merge(edges2,nodes2[,c(2,6,7,10,5)], by.x = "node_b", by.y = "node_id")
edges2$node_b<-NULL
edges2$node_a<-NULL
edges2$weight<-NULL


# nodes1<-nodes1[,c(11,6,7,10,5,9,8)]
# nodes1<-nodes1[,c(11,6,7,10,5,9,8)]
allNodes<-rbind(nodes1[,c(11,6,7,10,5,9,8)],nodes2[,c(11,6,7,10,5,9,8)])
allNodes<-allNodes[!duplicated(allNodes[,1:3]),]
qtdMapas=2
dimMapa=512
qtdNodes=max(nrow(nodes1),nrow(nodes2))
#qtdEdges=nrow(edges1)+nrow(edges2)
param<-"/home/clovis/projects/vcomplex/trunk/Dados/param.conf"
cat(file = param,sep = "\n",
    "#cor de fundo","1",
    "#quantidade de mapas",qtdMapas,
    "#dimensao do mapa",dimMapa,
    "#nodes",nrow(allNodes),
    "#edges",nrow(edges1),nrow(edges2))
allEdges<-rbind(edges1,edges2)
write.table(allNodes,
            file="/home/clovis/projects/vcomplex/trunk/Dados/nodes.csv",
            sep = ",",
            row.names = F,
            col.names = F,
            quote = F)
write.table(allEdges,
            file="/home/clovis/projects/vcomplex/trunk/Dados/edges.csv",
            sep = ",",
            row.names = F,
            col.names = F,
            quote = F)

#via complex

edgesVC<-edges[edges$node_a%in%nodes1$node_id &
                 edges$node_b%in%nodes1$node_id,]
edgesVC<-merge(edgesVC,nodes1[,c(2,6,7,11)], by.x = "node_a", by.y = "node_id")
edgesVC<-merge(edgesVC,nodes1[,c(2,6,7,11)], by.x = "node_b", by.y = "node_id")
edgesVC<-edgesVC[,c(6,9)]
edgesVC<-unique(edgesVC)


file="/home/clovis/.wine/drive_c/Arquivos de programas/ViaComplex/samples/netT1.dat"
cat(file = file, "*edges\t\n",append = F)
for(i in 1:nrow(edgesVC)){
  cat(file = file, 
      paste0(edgesVC$Symbol.x[i],"\t",
      edgesVC$Symbol.y[i],"\n"),
      append = T)
}
cat(file = file, "*nodes\t\t\n",append = T)

nodesVC<-nodes1
nodesVC<-nodesVC[,c(11,6,7)]
nodesVC$x<-(nodesVC$x-min(nodesVC$x))/(max(nodesVC$x)-min(nodesVC$x))
nodesVC$y<-(nodesVC$y-min(nodesVC$y))/(max(nodesVC$y)-min(nodesVC$y))
for(i in 1:nrow(nodesVC)){
  cat(file = file, 
      paste0(nodesVC$Symbol[i],"\t",
             nodesVC$x[i],"\t",
             nodesVC$y[i],"\n"),
      append = T)
}
file="/home/clovis/.wine/drive_c/Arquivos de programas/ViaComplex/samples/exp.dat"
cat(file = file, "ID\tGeneSymbol\tMCF7_High_met\tMCF7_Control\n",append = F)
expVC<-nodes1[,c(11,10)]
expVC$ct<-0
for(i in 1:nrow(expVC)){
  cat(file = file, 
      paste0(i,"\t",
             expVC$Symbol[i],"\t",
             expVC$logFC[i],"\t",
             expVC$ct[i],"\n"),
      append = T)
}
