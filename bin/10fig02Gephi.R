rm(list = ls())

lightColor = function(color){
  if(class(color) != "character"){
    cat("Color must be a string.")
  }
  R<-paste0("0x",substr(color,2,3))
  G<-paste0("0x",substr(color,4,5))
  B<-paste0("0x",substr(color,6,7))
  r<-ifelse(strtoi(R)<=127,strtoi(R)+128,255)
  g<-ifelse(strtoi(G)<=127,strtoi(G)+128,255)
  b<-ifelse(strtoi(B)<=127,strtoi(B)+128,255)
  color2<-paste0("#",
                sprintf("%x",r),
                sprintf("%x",g),
                sprintf("%x",b))
  return(color2)
}

setwd("/home/clovis/Dropbox/Chumbo/")
load("./transcriptograms/allTranscriptogramers80")
load("./Data/colors.RData")


gephiDir<-"/home/clovis/Dropbox/Chumbo/figuras/Figura4/redes/toGephi/"

setwd("/home/clovis/Dropbox/Chumbo/figuras/Figura4/redes/NodesEdges/")

# tempo 1 regra ----
#clusteres únicos
tempo<-1
ncl<-1
allCl<-c(1,2,3,4,5,7,9,10,11)
#string com os nomes corretos dos arquivos
nncl<-c("1_1","1_2","1_3","2","4",NA,"6",NA,"8_1","8_2","9")
for (ncl in allCl){
  nodes<-read.table(file = paste0("t",tempo,"c",ncl,"Nodes.txt"),
                    sep="\t",
                    header = T, 
                    stringsAsFactors = F)
  #normaliza coordenadas dos nodes entre 50 e 974
  nodes$x<-(nodes$x-min(nodes$x))/(max(nodes$x)-min(nodes$x))*974+50
  nodes$y<-(nodes$y-min(nodes$y))/(max(nodes$y)-min(nodes$y))*974+50
  
  nodes$Protein<-substr(nodes$alias,1,20)
  nodes$net<-substr(nodes$alias,24,25)
  nodes$alias<-NULL
  
  #parametros de cor
  nc0<-substr(toupper(color[[tempo]][ncl]),1,7)
  #cinza
  ncn<-"#AAAAAA"
  
  #acerta cor dos nós
  nodes$color[nodes$net==0]<-nc0

  #acerta parametros para o multigravity
  nodes$gravity_x[nodes$net==0]<-500
  nodes$gravity_y[nodes$net==0]<-500

  colnames(nodes)<-c("Id","shape","size","color","x","y","weight","Protein","net","gravity_x","gravity_y" )

  edges<-read.table(file = paste0("t",tempo,"c",ncl,"Edges.txt"),
                    sep="\t",
                    header = T, 
                    stringsAsFactors = F)
  edges<-merge(edges,nodes[,c(1,4,9)], by.x = "node_a", by.y = "Id")
  colnames(edges)<-c("Source","Target","weight","color","net")
  
  n0<-subset(x=nodes,nodes$net==0, select = Id)
  #check
  nrow(n0)
  #pesos de edge diferentes caso  pertençam a redes diferentes
  #rede pai será mais relaxada com weight = 5
  edges$weight<-0.001
  edges$weight[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-0.05

  #color = a de ambos os nós
  #se liga redes diferentes recebe cinza
  edges$color<-ncn
  edges$color[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-nc0

  edges<-unique(edges)
  write.csv(nodes,file=paste0(gephiDir,
                              "t",tempo,
                              "c",nncl[ncl],
                              "NodesGephi.csv"),
            row.names = F)
  write.csv(edges,file=paste0(gephiDir,
                              "t",tempo,
                              "c",nncl[ncl],
                              "EdgesGephi.csv"),
            row.names = F)
}  
  
# tempo 2 regra ----
#clusteres 1 pra 1
tempo<-2
allCl2<-c(2,4,6, 9)
allCl1<-c(4,5,7,11)
i=1
for (i in 1:length(allCl2)){
  ncl<-allCl2[i]
  ncl1<-allCl1[i]
  nodes<-read.table(file = paste0("t",tempo,"c",ncl,"Nodes.txt"),
                    sep="\t",
                    header = T, 
                    stringsAsFactors = F)
  #normaliza coordenadas dos nodes entre 50 e 974
  nodes$x<-(nodes$x-min(nodes$x))/(max(nodes$x)-min(nodes$x))*974+50
  nodes$y<-(nodes$y-min(nodes$y))/(max(nodes$y)-min(nodes$y))*974+50
  
  nodes$Protein<-substr(nodes$alias,1,20)
  nodes$net<-substr(nodes$alias,24,25)
  nodes$alias<-NULL
  
  #parametros de cor
  #teste
   # for(j in 1:11){
   #   #cat(j)
   #   barplot(c(seq(1:2)),
   #        col = c(toupper(color[[2]][j]),
   #                lightColor(toupper(color[[2]][j]))),
   #        names.arg = seq(1,2))
   # }
  nc0<-lightColor(toupper(color[[tempo]][allCl2[i]]))
  nc1<-substr(toupper(color[[1]][allCl1[i]]),1,7)
  barplot(c(seq(1:2)),
         col = c(nc0,nc1),
         names.arg = seq(1,2))

  #cinza
  ncn<-"#AAAAAA"
  
  #acerta cor dos nós
  nodes$color[nodes$net==0]<-nc0
  nodes$color[nodes$net==1]<-nc1
  #acerta parametros para o multigravity
  nodes$gravity_x[nodes$net==0]<-500
  nodes$gravity_y[nodes$net==0]<-500
  nodes$gravity_x[nodes$net==1]<-300
  nodes$gravity_y[nodes$net==1]<-250

  colnames(nodes)<-c("Id","shape","size","color","x","y","weight","Protein","net","gravity_x","gravity_y" )
  
  edges<-read.table(file = paste0("t",tempo,"c",ncl,"Edges.txt"),
                    sep="\t",
                    header = T, 
                    stringsAsFactors = F)
  edges<-merge(edges,nodes[,c(1,4,9)], by.x = "node_a", by.y = "Id")
  colnames(edges)<-c("Source","Target","weight","color","net")
  
  n0<-subset(x=nodes,nodes$net==0, select = Id)
  n1<-subset(x=nodes,nodes$net==1, select = Id)
  #check
  nrow(nodes)-nrow(n0)-nrow(n1)
  #pesos de edge diferentes caso  pertençam a redes diferentes
  #rede pai será mais relaxada com weight = 5
  edges$weight<-0.001
  edges$weight[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-0.05
  edges$weight[edges$Source%in%n1$Id & edges$Target%in%n1$Id]<-0.10
  #color = a de ambos os nós
  #se liga redes diferentes recebe cinza
  edges$color<-ncn
  edges$color[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-nc0
  edges$color[edges$Source%in%n1$Id & edges$Target%in%n1$Id]<-nc1

  edges<-unique(edges)

  edges<-unique(edges)
  edgesAux<-unique(c(edges$Source,edges$Target))
  
  nodesIn<-nodes[nodes$Id%in%edgesAux,]
  nodesOut<-nodes[!nodes$Id%in%edgesAux,]
  write.csv(nodesIn,file=paste0(gephiDir,            
                                "t",tempo,
                                "c",ncl,
                                "NodesInGephi.csv"),
            row.names = F)
  write.csv(nodesOut,file=paste0(gephiDir,            
                                 "t",tempo,
                                 "c",ncl,
                                 "NodesOutGephi.csv"),
            row.names = F)
  write.csv(edges,file=paste0(gephiDir,
                              "t",tempo,
                              "c",ncl,
                              "EdgesGephi.csv"),
            row.names = F)
}  

# tempo 2 cl 1 ----
tempo=2
ncl=1
nodes<-read.table(file = paste0("t",tempo,"c",ncl,"Nodes.txt"),
                  sep="\t",
                  header = T, 
                  stringsAsFactors = F)
#normaliza coordenadas dos nodes entre 50 e 974
nodes$x<-(nodes$x-min(nodes$x))/(max(nodes$x)-min(nodes$x))*974+50
nodes$y<-(nodes$y-min(nodes$y))/(max(nodes$y)-min(nodes$y))*974+50

nodes$Protein<-substr(nodes$alias,1,20)
nodes$net<-substr(nodes$alias,24,25)
nodes$alias<-NULL

#parametros de cor
nc0<-lightColor(toupper(color[[tempo]][ncl]))
nc1<-substr(toupper(color[[1]][1]),1,7)
nc2<-substr(toupper(color[[1]][2]),1,7)
nc3<-substr(toupper(color[[1]][3]),1,7)
#cinza
ncn<-"#AAAAAA"

#acerta cor dos nós
nodes$color[nodes$net==0]<-nc0
nodes$color[nodes$net==1]<-nc1
nodes$color[nodes$net==2]<-nc2
nodes$color[nodes$net==3]<-nc3
#acerta parametros para o multigravity
nodes$gravity_x[nodes$net==0]<-500
nodes$gravity_y[nodes$net==0]<-500
nodes$gravity_x[nodes$net==1]<-660
nodes$gravity_y[nodes$net==1]<-352
nodes$gravity_x[nodes$net==2]<-270
nodes$gravity_y[nodes$net==2]<-270
nodes$gravity_x[nodes$net==3]<-412
nodes$gravity_y[nodes$net==3]<-750

colnames(nodes)<-c("Id","shape","size","color","x","y","weight","Protein","net","gravity_x","gravity_y" )

edges<-read.table(file = paste0("t",tempo,"c",ncl,"Edges.txt"),
                  sep="\t",
                  header = T, 
                  stringsAsFactors = F)
edges<-merge(edges,nodes[,c(1,4,9)], by.x = "node_a", by.y = "Id")
colnames(edges)<-c("Source","Target","weight","color","net")

n0<-subset(x=nodes,nodes$net==0, select = Id)
n1<-subset(x=nodes,nodes$net==1, select = Id)
n2<-subset(x=nodes,nodes$net==2, select = Id)
n3<-subset(x=nodes,nodes$net==3, select = Id)
nrow(n0)+nrow(n1)+nrow(n2)+nrow(n3)
#pesos de edge diferentes caso  pertençam a redes diferentes
#rede pai será mais relaxada com weight = 5
edges$weight<-0.001
edges$weight[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-0.05
edges$weight[edges$Source%in%n1$Id & edges$Target%in%n1$Id]<-0.10
edges$weight[edges$Source%in%n2$Id & edges$Target%in%n2$Id]<-0.10
edges$weight[edges$Source%in%n3$Id & edges$Target%in%n3$Id]<-0.10
#color = a de ambos os nós
#se liga redes diferentes recebe cinza
edges$color<-ncn
edges$color[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-nc0
edges$color[edges$Source%in%n1$Id & edges$Target%in%n1$Id]<-nc1
edges$color[edges$Source%in%n2$Id & edges$Target%in%n2$Id]<-nc2
edges$color[edges$Source%in%n3$Id & edges$Target%in%n3$Id]<-nc3

edges<-unique(edges)
edgesAux<-unique(c(edges$Source,edges$Target))

nodesIn<-nodes[nodes$Id%in%edgesAux,]
nodesOut<-nodes[!nodes$Id%in%edgesAux,]
write.csv(nodesIn,file=paste0(gephiDir,"t2c1NodesInGephi.csv"),
          row.names = F)
write.csv(nodesOut,file=paste0(gephiDir,"t2c1NodesOutGephi.csv"),
          row.names = F)
write.csv(edges,file=paste0(gephiDir,"t2c1EdgesGephi.csv"),
          row.names = F)



#t2 só 3, 10 e 11
# tempo 2 únicos ----
#clusteres únicos
tempo<-2
ncl<-1
allCl<-c(3,5,7,10,11)
for (ncl in allCl){
  nodes<-read.table(file = paste0("t",tempo,"c",ncl,"Nodes.txt"),
                    sep="\t",
                    header = T, 
                    stringsAsFactors = F)
  #normaliza coordenadas dos nodes entre 50 e 974
  nodes$x<-(nodes$x-min(nodes$x))/(max(nodes$x)-min(nodes$x))*974+50
  nodes$y<-(nodes$y-min(nodes$y))/(max(nodes$y)-min(nodes$y))*974+50
  
  nodes$Protein<-substr(nodes$alias,1,20)
  nodes$net<-substr(nodes$alias,24,25)
  nodes$alias<-NULL
  
  #parametros de cor
  nc0<-substr(toupper(color[[tempo]][ncl]),1,7)
  #cinza
  ncn<-"#AAAAAA"
  
  #acerta cor dos nós
  nodes$color[nodes$net==0]<-nc0
  
  #acerta parametros para o multigravity
  nodes$gravity_x[nodes$net==0]<-500
  nodes$gravity_y[nodes$net==0]<-500
  
  colnames(nodes)<-c("Id","shape","size","color","x","y","weight","Protein","net","gravity_x","gravity_y" )
  
  edges<-read.table(file = paste0("t",tempo,"c",ncl,"Edges.txt"),
                    sep="\t",
                    header = T, 
                    stringsAsFactors = F)
  edges<-merge(edges,nodes[,c(1,4,9)], by.x = "node_a", by.y = "Id")
  colnames(edges)<-c("Source","Target","weight","color","net")
  
  n0<-subset(x=nodes,nodes$net==0, select = Id)
  #check
  nrow(n0)
  #pesos de edge diferentes caso  pertençam a redes diferentes
  #rede pai será mais relaxada com weight = 5
  edges$weight<-0.001
  edges$weight[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-0.05
  
  #color = a de ambos os nós
  #se liga redes diferentes recebe cinza
  edges$color<-ncn
  edges$color[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-nc0
  
  edges<-unique(edges)
  
  write.csv(nodes,file=paste0(gephiDir,
                              "t",tempo,
                              "c",ncl,
                              "NodesGephi.csv"),
            row.names = F)
  write.csv(edges,file=paste0(gephiDir,
                              "t",tempo,
                              "c",ncl,
                              "EdgesGephi.csv"),
            row.names = F)
}  

#t2 invertido t2c5 com t1c6, t2c7 com t1c8
# tempo 2 invertido ----
tempo<-2
allCl2<-c(5,7)
allCl1<-c(6,8)
i=1
for (i in 1:length(allCl2)){
  ncl<-allCl2[i]
  ncl1<-allCl1[i]
  nodes<-read.table(file = paste0("t",tempo,"c",ncl,"Nodes.txt"),
                    sep="\t",
                    header = T, 
                    stringsAsFactors = F)
  #normaliza coordenadas dos nodes entre 50 e 974
  nodes$x<-(nodes$x-min(nodes$x))/(max(nodes$x)-min(nodes$x))*974+50
  nodes$y<-(nodes$y-min(nodes$y))/(max(nodes$y)-min(nodes$y))*974+50
  
  nodes$Protein<-substr(nodes$alias,1,20)
  nodes$net<-substr(nodes$alias,24,25)
  nodes$alias<-NULL
  
  #parametros de cor
  #teste
  # for(j in 1:11){
  #   #cat(j)
  #   barplot(c(seq(1:2)),
  #        col = c(toupper(color[[2]][j]),
  #                lightColor(toupper(color[[2]][j]))),
  #        names.arg = seq(1,2))
  # }
  nc0<-lightColor(toupper(color[[tempo]][allCl2[i]]))
  nc1<-substr(toupper(color[[1]][allCl1[i]]),1,7)
  barplot(c(seq(1:2)),
          col = c(nc0,nc1),
          names.arg = seq(1,2))
  
  #cinza
  ncn<-"#AAAAAA"
  
  #acerta cor dos nós
  nodes$color[nodes$net==0]<-nc0
  nodes$color[nodes$net==1]<-nc1
  #acerta parametros para o multigravity
  nodes$gravity_x[nodes$net==1]<-500
  nodes$gravity_y[nodes$net==1]<-500
  nodes$gravity_x[nodes$net==0]<-300
  nodes$gravity_y[nodes$net==0]<-250
  
  colnames(nodes)<-c("Id","shape","size","color","x","y","weight","Protein","net","gravity_x","gravity_y" )
  
  edges<-read.table(file = paste0("t",tempo,"c",ncl,"Edges.txt"),
                    sep="\t",
                    header = T, 
                    stringsAsFactors = F)
  edges<-merge(edges,nodes[,c(1,4,9)], by.x = "node_a", by.y = "Id")
  colnames(edges)<-c("Source","Target","weight","color","net")
  
  n0<-subset(x=nodes,nodes$net==0, select = Id)
  n1<-subset(x=nodes,nodes$net==1, select = Id)
  #check
  nrow(nodes)-nrow(n0)-nrow(n1)
  #pesos de edge diferentes caso  pertençam a redes diferentes
  #rede pai será mais relaxada com weight = 5
  edges$weight<-0.001
  edges$weight[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-0.05
  edges$weight[edges$Source%in%n1$Id & edges$Target%in%n1$Id]<-0.10
  #color = a de ambos os nós
  #se liga redes diferentes recebe cinza
  edges$color<-ncn
  edges$color[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-nc0
  edges$color[edges$Source%in%n1$Id & edges$Target%in%n1$Id]<-nc1
  
  edges<-unique(edges)

  edgesAux<-unique(c(edges$Source,edges$Target))
  
  nodesIn<-nodes[nodes$Id%in%edgesAux,]
  nodesOut<-nodes[!nodes$Id%in%edgesAux,]
  write.csv(nodesIn,file=paste0(gephiDir,            
                                "t",tempo,
                                "c",ncl,
                                "NodesInGephi.csv"),
            row.names = F)
  write.csv(nodesOut,file=paste0(gephiDir,            
                                 "t",tempo,
                                 "c",ncl,
                                 "NodesOutGephi.csv"),
            row.names = F)
  write.csv(edges,file=paste0(gephiDir,
                              "t",tempo,
                              "c",ncl,
                              "EdgesGephi.csv"),
            row.names = F)
  
}  

#t2 contem t2c8 contem t2c9 e t2c10 
# tempo 2 cl 8 ----
tempo=2
ncl=8
nodes<-read.table(file = paste0("t",tempo,"c",ncl,"Nodes.txt"),
                  sep="\t",
                  header = T, 
                  stringsAsFactors = F)
#normaliza coordenadas dos nodes entre 50 e 974
nodes$x<-(nodes$x-min(nodes$x))/(max(nodes$x)-min(nodes$x))*974+50
nodes$y<-(nodes$y-min(nodes$y))/(max(nodes$y)-min(nodes$y))*974+50

nodes$Protein<-substr(nodes$alias,1,20)
nodes$net<-substr(nodes$alias,24,25)
nodes$alias<-NULL

#parametros de cor
nc0<-lightColor(toupper(color[[tempo]][ncl]))
nc1<-substr(toupper(color[[1]][9]),1,7)
nc2<-substr(toupper(color[[1]][10]),1,7)
#cinza
ncn<-"#AAAAAA"

#acerta cor dos nós
nodes$color[nodes$net==0]<-nc0
nodes$color[nodes$net==1]<-nc1
nodes$color[nodes$net==2]<-nc2
#acerta parametros para o multigravity
nodes$gravity_x[nodes$net==0]<-500
nodes$gravity_y[nodes$net==0]<-500
nodes$gravity_x[nodes$net==1]<-750
nodes$gravity_y[nodes$net==1]<-750
nodes$gravity_x[nodes$net==2]<-250
nodes$gravity_y[nodes$net==2]<-250

colnames(nodes)<-c("Id","shape","size","color","x","y","weight","Protein","net","gravity_x","gravity_y" )

edges<-read.table(file = paste0("t",tempo,"c",ncl,"Edges.txt"),
                  sep="\t",
                  header = T, 
                  stringsAsFactors = F)
edges<-merge(edges,nodes[,c(1,4,9)], by.x = "node_a", by.y = "Id")
colnames(edges)<-c("Source","Target","weight","color","net")

n0<-subset(x=nodes,nodes$net==0, select = Id)
n1<-subset(x=nodes,nodes$net==1, select = Id)
n2<-subset(x=nodes,nodes$net==2, select = Id)
nrow(n0)+nrow(n1)+nrow(n2)
#pesos de edge diferentes caso  pertençam a redes diferentes
#rede pai será mais relaxada com weight = 5
edges$weight<-0.001
edges$weight[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-0.05
edges$weight[edges$Source%in%n1$Id & edges$Target%in%n1$Id]<-0.10
edges$weight[edges$Source%in%n2$Id & edges$Target%in%n2$Id]<-0.10
#color = a de ambos os nós
#se liga redes diferentes recebe cinza
edges$color<-ncn
edges$color[edges$Source%in%n0$Id & edges$Target%in%n0$Id]<-nc0
edges$color[edges$Source%in%n1$Id & edges$Target%in%n1$Id]<-nc1
edges$color[edges$Source%in%n2$Id & edges$Target%in%n2$Id]<-nc2

edges<-unique(edges)
edgesAux<-unique(c(edges$Source,edges$Target))

nodesIn<-nodes[nodes$Id%in%edgesAux,]
nodesOut<-nodes[!nodes$Id%in%edgesAux,]
write.csv(nodesIn,file=paste0(gephiDir,            
                              "t",tempo,
                              "c",ncl,
                              "NodesInGephi.csv"),
          row.names = F)
write.csv(nodesOut,file=paste0(gephiDir,            
                               "t",tempo,
                               "c",ncl,
                               "NodesOutGephi.csv"),
          row.names = F)
write.csv(edges,file=paste0(gephiDir,
                            "t",tempo,
                            "c",ncl,
                            "EdgesGephi.csv"),
          row.names = F)
