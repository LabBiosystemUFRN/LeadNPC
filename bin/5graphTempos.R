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
#load(file = "./Data/clusteres.RData")
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

for(i in 1:max(c1)){
    clusters12$Clust1[clusters12$Position%in%c(clusters[[1]]$ini[c1[i]]:clusters[[1]]$fim[c1[i]])]<-c1[i]
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
save(cl12,file = "./Data/cl12SB.RData")
load(file = "./Data/cl12SB.RData")
#Grapho master ----
#cria grafo base para os outros

gMaster <- graph_from_data_frame(transc[[1]]@association, directed = F)
gMaster<-att.mapv(g = gMaster, dat = cl12, refcol = 1)
gMaster<-att.setv(g = gMaster, from = "Symbol", to = "nodeAlias" )
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
  #inicioLoop ----
  ref2=4
  myTheme$nestLineType="DOTTED"
  
  for(ref2 in 1:max(c2)){
    #SetRef----  
    #ref2=13 
    ref1<-cltRef$c1[cltRef$c2 ==  ref2]
  ref1<-na.exclude(ref1)
  #ref5<-unique(cltRef$c5[cltRef$c4%in%ref4])
  
  genes1<-cl12[cl12$Clust1%in%ref1,]
  genes2<-cl12[cl12$Clust2%in%ref2,]

  sgenes0<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$nodeAlias %in% genes1$Protein])
  
  #Intervalo1 ----
  sgenes1<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% genes1$Protein])
  if(ref2 == 3){
    V(sgenes1)$nodeLineColor<-"#999999"
  }else{
    V(sgenes1)$nodeLineColor<-color[[2]][ref2]
    }
  V(sgenes1)$nodeColor <- color[[2]][ref2]
  V(sgenes1)$size <- 10
  name1<-paste0("Cluster",ref2,"Intervalo 1")
  
  resetd(rdp)
  
  myTheme$nestAlias<-name1
  myTheme$nestLineColor<-"#bfbfefFF"
  
  RedeR::addGraph(rdp, sgenes1,
                  isNest=T,
                  gcoord=c(50,50), gscale=200, theme=myTheme,zoom=60)
  deSelectNodes(rdp)
  relax(rdp,70,10,150,150,150,75,5,50,ps=T)
  
  readline(prompt=paste0("C",ref2, "T1. Enter continua... "))
  
  
  #Intervalo2 ----
  sgenes2<-subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% genes2$Protein])
  if(ref2 == 3){
    V(sgenes2)$nodeLineColor<-"#999999"
  }else{
    V(sgenes2)$nodeLineColor<-color[[2]][ref2]
  }
  #V(sgenes2)$nodeLineColor<-color[[2]][ref2]
  # V(sgenes2)$nodeColor <- ifelse(V(sgenes2)$nodeAlias %in% tfs, 
  #                                "#666666FF", #Fator de transcriçaõ
  #                                color[[2]][ref2]) #cor 2
  # V(sgenes2)$size <- ifelse(V(sgenes2)$nodeAlias %in% tfs, 20, 10)
  
  V(sgenes2)$nodeColor <- color[[2]][ref2] #cor 2
  V(sgenes2)$size <- 10
  
    sgenes2s<-subgraph(graph=sgenes2, v=V(sgenes2)[V(sgenes2)$name %in% genes1$Protein])
    name2<-paste0("Cluster",ref2,"Intervalo 2")

    escala<-round(sum(genes2$Protein%in%genes1$Protein)/nrow(genes2)*100,0)
    resetd(rdp)
    myTheme$nestAlias<-name2
    myTheme$nestLineColor<-"#bfbfefFF"
    RedeR::addGraph(rdp, sgenes2,
                    isNest=T,
                    gcoord=c(50,50), gscale=200, theme=myTheme,zoom=40)
    myTheme$nestAlias<-""
    myTheme$nestLineColor<-"#444444FF"
    nestNodes(rdp, nodes = V(sgenes2s)$name, 
              #gcoord=c(25,25),
              parent="N0", 
              gscale=escala,
              theme=myTheme)

    deSelectNodes(rdp)
    relax(rdp,70,10,150,150,150,75,5,50,ps=T)

    readline(prompt=paste0("C",ref2, "T2. Enter continua... "))
    
}
#fim############################################################################  
    rdp <- RedPort()
    RedeR::calld(rdp)
    
    assoc<-transc[[3]]@association
    assoc<-assoc[assoc$p1%in%cl12$Protein & assoc$p2%in%cl12$Protein,]
    gMaster <- graph_from_data_frame(assoc, directed = F)
    
    RedeR::addGraph(rdp, gMaster,
                    isNest=F,
                    gcoord=c(50,50), gscale=100, theme=myTheme,zoom=10)
    deSelectNodes(rdp)
    relax(rdp,70,10,100,150,150,75,5,20,ps=F)
    
    save(gMaster,file = "gMaster.RDat")

    
    