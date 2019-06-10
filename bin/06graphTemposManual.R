rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadNPC/"
setwd(baseDir)
source(file = "bin/00base.R")

load("./Data/allTranscriptogramers80")
load("./Data/tfs.RData")
load("./Data/colors.RData")

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



#intersecção dos clusteres
 load(file = "./Data/clusteres.RData")
#load(file = "./Data/clusteresSB.RData")
c1<-clusteres$T1
c2<-clusteres$T2

clusters12<-data.frame(Protein = transc[[2]]@ordering$Protein,
                        Position = transc[[2]]@ordering$Position,
                        Clust1 = NA,
                        Clust2 = NA)

for(i in 1:12){
    clusters12$Clust1[clusters12$Position%in%c(clusters[[1]]$ini[i]:clusters[[1]]$fim[i])]<-i
}
for(i in 1:12){
  clusters12$Clust2[clusters12$Position%in%c(clusters[[2]]$ini[i]:clusters[[2]]$fim[i])]<-i
}
clusters12$Clust1[clusters12$Clust1==12]<-1
clusters12$Clust2[clusters12$Clust2==12]<-1


# save(clusters12,file='./Data/clusters12.RData')
 load(file='./Data/clusters12.RData')

 protSymbol<-transc[[1]]@Protein2Symbol

colnames(protSymbol)<-c("Protein", "Symbol")
clusters12<-merge(clusters12,protSymbol, by="Protein", all.x = T)
clusters12$tFactor<-"ELLIPSE"
#clusters12$tFactor[clusters12$Symbol%in%tfs]<-"DIAMOND"
clusters12$Cor1<-color[[1]][clusters12$Clust1]
clusters12$Cor2<-color[[2]][clusters12$Clust2]

cl12<-(clusters12[!(is.na(clusters12$Clust1) & is.na(clusters12$Clust2)),])
 load(file = "./Data/cl12.RData")
 load(file='./Data/clusters12.RData')
#Grapho master ----
#cria grafo base para os outros

gMaster <- graph_from_data_frame(transc[[1]]@association, directed = F)
gMaster<-att.mapv(g = gMaster, dat = cl12, refcol = 1)
V(gMaster)$nodeShape <- "CIRCLE"
E(gMaster)$edgeColor<-"#c9c9c9ff"

#cria dtf com intersecção
cltRef<-data.frame(c1=c1,c2=c2)

#inicia o port
rdp <- RedPort()
RedeR::calld(rdp)
  #extrai genes para grupos 3 4 5 com base no 3 ref3
  #inicioLoop ----
  ref2=4
  myTheme$nestLineType="DOTTED"

  # Cluster 1 int 2 ----
  #rede pai
  nodesP<-subset(clusters12,Clust2==1,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust1==1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  
  nodesF2<-subset(clusters12,Clust1==2,select = c("Protein") )
  nodesF2$Protein<-as.character(nodesF2$Protein)

  nodesF3<-subset(clusters12,Clust1==3,select = c("Protein") )
  nodesF3$Protein<-as.character(nodesF3$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[2]][1]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[1]][1]
  V(pai)$nodeColor[V(pai)$name %in% nodesF2$Protein]<-color[[1]][2]
  V(pai)$nodeColor[V(pai)$name %in% nodesF3$Protein]<-color[[1]][3]
  V(pai)$size <- 10
  #name1<-paste0("Cluster",ref2,"Intervalo 1")



  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#bfbfefFF"
  resetd(rdp)
  RedeR::addGraph(rdp, pai,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=30)
  myTheme$nestAlias<-""    

  
  myTheme$nestLineColor<-"#bfbfefFF"
  myTheme$nestAlias<-""
  myTheme$nestLineColor<-"#444444FF"
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*200*2,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  escala<-round(sum(nodesP$Protein%in%nodesF2$Protein)/nrow(nodesP)*200*2,0)
  nestNodes(rdp, nodes = nodesF2$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  escala<-round(sum(nodesP$Protein%in%nodesF3$Protein)/nrow(nodesP)*200*2,0)
  nestNodes(rdp, nodes = nodesF3$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  
  deSelectNodes(rdp)
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  
  # Cluster 2 int 2 ----
  #rede pai
  nPai=2
  nF1=4
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
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  

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
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  
  # Cluster 4 int 2 ----
  #rede pai
  nPai=4
  nF1=5
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
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  
  # Cluster 6 int 2 ----
  #rede pai
  nPai=6
  nF1=7
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
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*100,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)

  deSelectNodes(rdp)
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  
  # Cluster 7 int 2 ----
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
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  
  # Cluster 8 int 2 ----
  #rede pai
  nPai=8
  nF1=9
  nF2=10
  nodesP<-subset(clusters12,Clust2==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust1==nF1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  nodesF2<-subset(clusters12,Clust1==nF2,select = c("Protein") )
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
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  RedeR::selectNodes(rdp,nodesF1$Protein)
  
  # Cluster 9 int 2 ----
  #rede pai
  nPai=9
  nF1=11
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
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
  
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

  deSelectNodes(rdp)
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)

  
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
  
  deSelectNodes(rdp)
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  
  
  # Cluster 5 int 2 ----
  #está invertido. cluster 6 é maior que o 5
  #rede pai
  nPai=6
  nF1=5
  nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
  nodesP$Protein<-as.character(nodesP$Protein)
  
  nodesF1<-subset(clusters12,Clust2==nF1,select = c("Protein") )
  nodesF1$Protein<-as.character(nodesF1$Protein)
  
  pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
  
  #cores
  V(pai)$nodeColor<-color[[1]][nPai]
  V(pai)$nodeLineColor<-"#999999"
  V(pai)$nodeColor[V(pai)$name %in% nodesF1$Protein]<-color[[2]][nF1]
  V(pai)$size <- 10

  
  
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
  escala<-round(sum(nodesP$Protein%in%nodesF1$Protein)/nrow(nodesP)*100,0)
  nestNodes(rdp, nodes = nodesF1$Protein, 
            #gcoord=c(25,25),
            parent="N0", 
            gscale=escala,
            theme=myTheme)
  
  deSelectNodes(rdp)
  #relax(rdp,70,10,150,150,150,75,150,50,ps=T)
  RedeR::selectNodes(rdp,nodesF1$Protein)
  
  # Todos Intervalo 1 ----
  nPai=7
  for(nPai in 1:11){
    nodesP<-subset(clusters12,Clust1==nPai,select = c("Protein") )
    nodesP$Protein<-as.character(nodesP$Protein)
    pai<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
    
    #cores
    V(pai)$nodeColor<-color[[1]][nPai]
    V(pai)$nodeLineColor<-"#555555"
    V(pai)$size <- 30

    
    myTheme$nestAlias<-""
    myTheme$nestLineColor<-"#444444FF"
    resetd(rdp)
    RedeR::addGraph(rdp, pai,
                    isNest=T,
                    gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
    myTheme$nestAlias<-""    
    
    deSelectNodes(rdp)
    #relax(rdp,70,10,150,150,150,75,150,50,ps=T)
    
    #aguarda
    readline(prompt=paste0("i1c",nPai, ". Enter continua... "))
    
  }
    
