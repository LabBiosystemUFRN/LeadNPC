rm(list = ls())
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

myTheme<-list()
zoom=30
myTheme$nestShape='ROUNDED_RECTANGLE'
myTheme$nestColor='#ffffff'
myTheme$nestLineWidth=4*(100/zoom)
myTheme$nestFontSize=24*(100/zoom)
myTheme$nestFontX=5
myTheme$nestFontY=10.8
myTheme$isAssign=TRUE


p<-matrix(0,ncol = 2, nrow = nrow(transc[[1]]@ordering))
i=1
for(i in 1:2){
  object<-transc[[i]]
  j=1
  for(j in 1:length(object@clusters)){
    pontos<-object@clusters[[j]]
    #colnames(pontos)<-c("x","y")
    p[pontos[[1]]:pontos[[2]],i]<-2^(i)
  }
}
genesIndex<-p[,1]+p[,2]
#se 1 G3, se 2 G4, se 4 G5, se 3 G3 e G4, 6 G4 e G5, se 7 G3 G4 e G5 se 0 ninguem

g0<-object@ordering$Protein[genesIndex==0]
g1<-object@ordering$Protein[genesIndex==2]
g12<-object@ordering$Protein[genesIndex==6]
g2<-object@ordering$Protein[genesIndex==4]

#check
length(g1)+length(g12)+length(g2)+length(g0)-
  nrow(object@ordering)

#limites dos clusteres
clusters<-list()
object<-transc[[1]]
for(i in 1:2){
  c1<-(unlist(map(transc[[i]]@clusters, 1)))
  c2<-(unlist(map(transc[[i]]@clusters, 2)))
  clusters[[i]]<-data.frame(ini=c1,fim=c2)
}
i=1

# #figura de referencia das intersecções
# p1<-ggplot()+
#   theme_bw()+
#   scale_y_continuous(limits = c(0,20))+
#   scale_x_continuous(limits = c(0,nrow(object@ordering)))+
#   xlab("Position")
# i=4
# for(i in 5:5){
#   dtf<-clusters[[i]]
#   j=5
#   for(j in 1:3){#nrow(clusters[[i]])){
#     p1 <- p1 + geom_segment(aes(x = clusters[[i]]$ini[j], 
#                                 y = i, 
#                                 xend = clusters[[i]]$fim[j], 
#                                 yend = i), col=(i-1),size=4)
#     summary(p1)
#   }
# }
# print(p1)




#intersecção dos clusteres
# load(file = "./Data/clusteres.RData")
load(file = "./Data/clusteresSB.RData")
c1<-clusteres$T1
c2<-clusteres$T2

clusters12<-data.frame(Protein = transc[[2]]@ordering$Protein,
                        Position = transc[[2]]@ordering$Position,
                        Clust1 = NA,
                        Clust2 = NA)
# i=5
# for(i in 1:length(c1)){
#   if(!is.na(c1[i])){
#     clusters12$Clust1[clusters12$Position%in%c(clusters[[1]]$ini[c1[i]]:clusters[[1]]$fim[c1[i]])]<-c1[i]
#   }
#   if(!is.na(c2[i])){
#     clusters12$Clust2[clusters12$Position%in%c(clusters[[2]]$ini[c2[i]]:clusters[[2]]$fim[c2[i]])]<-c2[i]
#   }
# }

############################3
#sem bound
# for(i in 1:12){
#     clusters12$Clust1[clusters12$Position%in%c(clusters[[1]]$ini[i]:clusters[[1]]$fim[i])]<-i
# }
# for(i in 1:12){
#   clusters12$Clust2[clusters12$Position%in%c(clusters[[2]]$ini[i]:clusters[[2]]$fim[i])]<-i
# }
# clusters12$Clust1[clusters12$Clust1==12]<-1
# clusters12$Clust2[clusters12$Clust2==12]<-1

for(i in 1:15){
  clusters12$Clust1[clusters12$Position%in%c(clusters[[1]]$ini[i]:clusters[[1]]$fim[i])]<-i
}
for(i in 1:26){
  clusters12$Clust2[clusters12$Position%in%c(clusters[[2]]$ini[i]:clusters[[2]]$fim[i])]<-i
}


#   if(!is.na(c2[i])){
#     clusters12$Clust2[clusters12$Position%in%c(clusters[[2]]$ini[c2[i]]:clusters[[2]]$fim[c2[i]])]<-c2[i]
#   }
# }


# save(clusters12,file='./Data/clusters12.RData')
# load(file='./Data/clusters12.RData')
save(clusters12,file='./Data/clusters12SB.RData')
load(file='./Data/clusters12SB.RData')
protSymbol<-transc[[1]]@Protein2Symbol

colnames(protSymbol)<-c("Protein", "Symbol")
clusters12<-merge(clusters12,protSymbol, by="Protein", all.x = T)
clusters12$tFactor<-"ELLIPSE"
#clusters12$tFactor[clusters12$Symbol%in%tfs]<-"DIAMOND"
clusters12$Cor1<-color[[1]][clusters12$Clust1]
clusters12$Cor2<-color[[2]][clusters12$Clust2]

cl12<-(clusters12[!(is.na(clusters12$Clust1) & is.na(clusters12$Clust2)),])
# save(cl12,file = "./Data/cl12.RData")
# load(file = "./Data/cl12.RData")
# load(file='./Data/clusters12.RData')
save(cl12,file = "./Data/cl12SB.RData")
load(file = "./Data/cl12SB.RData")
load(file='./Data/clusters12SB.RData')
#Grapho master ----
#cria grafo base para os outros

gMaster <- graph_from_data_frame(transc[[1]]@association, directed = F)
gMaster<-att.mapv(g = gMaster, dat = cl12, refcol = 1)
#gMaster<-att.setv(g = gMaster, from = "Symbol", to = "nodeAlias" )
#V(gMaster)$nodeShape <- ifelse(V(gMaster)$nodeAlias %in% tfs, "DIAMOND", "CIRCLE")
V(gMaster)$nodeShape <- "CIRCLE"
#V(gMaster)$nodeShape<-""
#gMaster<-att.setv(g = gMaster, from = "tFactor", to = "nodeShape" )
E(gMaster)$edgeColor<-"#c9c9c9ff"
#checagem do fator de transcrição
#check<-data.frame(sh=vertex_attr(gMaster,"nodeShape"),al=vertex_attr(gMaster,"nodeAlias"))
#check<-check[check$sh=="DIAMOND",]
#check[!(check$al%in%tfs),]
#nrow(check[(check$al%in%tfs),])-nrow(check)
#rm(check)

#intersecção dos clusteres
# c3<-c(1,2,3, 4, 5, 6, 7, 8)
# c4<-c(3,5,9,10,12,13,18,19)
# c5<-c(2,3,5, 6, 8,0,16,17)
#cria dtf com intersecção
cltRef<-data.frame(c1=c1,c2=c2)

#inicia o port
rdp <- RedPort()
RedeR::calld(rdp)
  #extrai genes para grupos 3 4 5 com base no 3 ref3
  ref2=4
  myTheme$nestLineType="DOTTED"


  #   #Intervalo2 ----
  # Cluster 1 int 2 ----
  #rede pai
  nPai=1
  nF1=1
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust1==nF1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[1]][nF1]
  V(pai)$size <- 10
  #name1<-paste0("Cluster",ref2,"Intervalo 1")
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*200,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  # Cluster 2 int 2 ----
  #rede pai
  nPai=2
  nF1=2
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust1==nF1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[1]][nF1]
  V(pai)$size <- 10
  #name1<-paste0("Cluster",ref2,"Intervalo 1")
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*200,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  

  # Cluster 3 int 2 ----
  #rede pai
  nPai=3
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)

  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10

  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"

  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  
  # Cluster 4 int 2 ----
  #rede pai
  nPai=4
  nF1=3
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust1==nF1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[1]][nF1]
  V(pai)$size <- 10

  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*200,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)

  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  # Cluster 5 int 2 ----
  #rede pai
  nPai=5
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 6 int 2 ----
  #rede pai
  nPai=6
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 7 int 2 ----
  #rede pai
  nPai=7
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 8 int 2 ----
  #rede pai
  nPai=8
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 9 int 2 ----
  #rede pai
  nPai=9
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  

  # Cluster 10 int 2 ----
  #rede pai
  nPai=10
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  
  # Cluster 11 int 2 ----
  #rede pai
  nPai=11
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 12 int 2 ----
  #rede pai
  nPai=12
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 13 int 2 ----
  #rede pai
  nPai=13
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 10 int 2 ----
  #rede pai
  nPai=14
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 15 int 2 ----
  #rede pai
  nPai=15
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 16 int 2 ----
  #rede pai
  nPai=16
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 17 int 2 ----
  #rede pai
  nPai=17
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 18 int 2 ----
  #rede pai
  nPai=18
  nF1=9
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust1==nF1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[1]][nF1]
  V(pai)$size <- 10
  #name1<-paste0("Cluster",ref2,"Intervalo 1")
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*200,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  
  # Cluster 19 int 2 ----
  #rede pai
  nPai=19
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 20 int 2 ----
  #rede pai
  nodesP<-subset(clusters12,Clust2==20,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)

  nodesF1<-subset(clusters12,Clust1==11,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)

  nodesF2<-subset(clusters12,Clust1==12,select = c("Protein") )
  nodesF2$Protein<-as.character(nodesF2$Protein)

  nodesF3<-subset(clusters12,Clust1==13,select = c("Protein") )
  nodesF3$Protein<-as.character(nodesF3$Protein)

  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])

  #cores
  V(pai)$nodeColor<-color[[2]][1]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[1]][11]
  V(pai)$nodeColor[V(pai)$name %in% nodesF2$Protein]<-color[[1]][12]
  V(pai)$nodeColor[V(pai)$name %in% nodesF3$Protein]<-color[[1]][13]
  V(pai)$size <- 10
  #name1<-paste0("Cluster",ref2,"Intervalo 1")



  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""


  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*200,0)
  nestNodes(rdp, nodes = nodesF1$Protein,
            #gcoord=c(25,25),
            parent="N0",
            gscale=escala,
            theme=myTheme)
  escala<-round(sum(nodesP$Protein%in%nodesF2$Protein)/nrow(nodesP)*200,0)
  nestNodes(rdp, nodes = nodesF2$Protein,
            #gcoord=c(25,25),
            parent="N0",
            gscale=escala,
            theme=myTheme)
  escala<-round(sum(nodesP$Protein%in%nodesF3$Protein)/nrow(nodesP)*200,0)
  nestNodes(rdp, nodes = nodesF3$Protein,
            #gcoord=c(25,25),
            parent="N0",
            gscale=escala,
            theme=myTheme)

  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  
  # Cluster 21 int 2 ----
  #rede pai
  nPai=21
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 22 int 2 ----
  #rede pai
  nPai=22
  nF1=14
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust1==nF1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[1]][nF1]
  V(pai)$size <- 10
  #name1<-paste0("Cluster",ref2,"Intervalo 1")
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*200,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  
  # Cluster 23 int 2 ----
  #rede pai
  nPai=23
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 24 int 2 ----
  #rede pai
  nPai=24
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 26 int 2 ----
  #rede pai
  nPai=26
  nF1=15
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust1==nF1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[1]][nF1]
  V(pai)$size <- 10
  #name1<-paste0("Cluster",ref2,"Intervalo 1")
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*200,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  
  #   #Intervalo1 ----
  # Cluster 1 int 1 ----
  #rede pai
  nPai=1
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  # Cluster 2 int 1 ----
  #rede pai
  nPai=2
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  # Cluster 3 int 1 ----
  #rede pai
  nPai=3
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  # Cluster 4 int 1 ----
  #rede pai
  nPai=4
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  # Cluster 5 int 1 ----
  #rede pai
  nPai=5
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  # Cluster 6 int 1 ----
  #rede pai
  nPai=6
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  
  # Cluster 7 int 1 ----
  #estaq invertido -  rede maior t1
  #rede pai
  nPai=7
  nF1=14
  nF2=15
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust2==nF1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  nodesF2<-subset(clusters12,Clust2==nF2,select = c("Protein") )
  nodesF2$Protein<-as.character(nodesF2$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[1]][nF1]
  V(pai)$nodeColor[V(pai)$name %in% nodesF2$Protein]<-color[[1]][nF2]
  V(pai)$size <- 10
  #name1<-paste0("Cluster",ref2,"Intervalo 1")
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*100,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  escala<-round(sum(nodesP$Protein%in%nodesF2$Protein)/nrow(nodesP)*100,0)
  nestNodes(rdp, nodes = nodesF2$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 8 int 2 ----
  #está invertido. cluster 8 é maior que o 7
  #rede pai
  nPai=8
  nF1=7
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust2==nF1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][7]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[2]][7]
  V(pai)$size <- 10
  #name1<-paste0("Cluster",ref2,"Intervalo 1")
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  escala<-85#round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*100,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 9 int 1 ----
  #rede pai
  nPai=9
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 10 int 1 ----
  #rede pai
  nPai=10
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 11 int 1 ----
  #rede pai
  nPai=11
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 12 int 1 ----
  #rede pai
  nPai=12
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 13 int 1 ----
  #rede pai
  nPai=13
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 14 int 1 ----
  #rede pai
  nPai=14
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 15 int 1 ----
  #rede pai
  nPai=15
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$size <- 10
  
  
  
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
  myTheme$nestAlias<-""    
  
  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
