rm(list = ls())

library("purrr")
library(ggplot2)

setwd("/home/clovis/Dropbox/Chumbo/")
#load("./Data/pheno_data.RData")
load("./Data/counts.RData")
load("./transcriptograms/allTranscriptogramers80")

figuras="figuras"

clusters<-list()
object<-transc[[1]]
for(i in 1:2){
  c1<-(unlist(map(transc[[i]]@clusters, 1)))
  c2<-(unlist(map(transc[[i]]@clusters, 2)))
  clusters[[i]]<-data.frame(ini=c1,fim=c2)
}

#figura de referencia das intersecções  ----
p1<-ggplot()+
  theme_bw()+
  scale_y_continuous(limits = c(0,20))+
  scale_x_continuous(limits = c(0,nrow(object@ordering)))+
  xlab("Position")
i=4

load("./Data/colors.RData")
for(i in 1:2){
  dtf<-clusters[[i]]
  j=5
  for(j in 1:nrow(clusters[[i]])){
    x1=dtf[j,1]
    x2=dtf[j,2]
    y1=i*2
    y2=i*2
    p1 <- p1 + geom_line(aes_string(x = c(x1,x2),y=c(y1,y2)), col=color[[i]][j],size=20)+
      annotate(geom="text", x=(x1+(x2-x1)/2), y=i*1.5, 
               label=j,
               color="black",cex=2)
    #summary(p1)
  }
}
svg(file = paste0("./",figuras,"/superposClust2.svg"))
plot(p1)
dev.off()


i=1
numLin<-nrow(transc[[1]]@ordering)
fn <- "/home/clovis/Dropbox/Chumbo/clusters.txt"
if (file.exists(fn)) file.remove(fn)

#cria lista com limites dos clusteres
pBreaks<-list()
#next cluster a ser processado
nxt<-c(0,0)
#primeiro cluster é circular?
crc<-c(0,0)
#cluster máximo
top<-c(0,0)
for(i in 1:2){
  object <- transc[[i]]
  nClust<-length(object@clusters)
  tmp<-data.frame(clNo=seq(1:nClust),
                  ini=sapply(object@clusters, function(x){x[1]}),
                  fim=sapply(object@clusters, function(x){x[2]}))
  pBreaks[[i]]<- tmp
  #inicio do processaemnto dos clusteres
  #0 indica que ultimo cluster continua no primeiro
  if(tmp$ini[1] == 0){
    nxt[i]<-1
    crc[i]<-1
    top[i]<-nClust
  }else{
    nxt[i]<-1
    crc[i]<-0
    top[i]<-nClust
  }
  rm(tmp,nClust)
}

#código para resultados
# 0-não superpostos
# 1-contido
# 2-contém
# 3-superposição head
# 4-superposição tail
#numero depois da virgula indica percentual de sobreposição
int1to2<-data.frame(cl1=numeric(),
                cl2=numeric(), 
                rel=numeric(),
                perc=numeric())
int2to1<-data.frame(cl1=numeric(),
                  cl2=numeric(), 
                  rel=numeric(),
                  perc=numeric())

correspond<-data.frame(cl1=numeric(),
                    cl2=numeric(), 
                    rel=numeric(),
                    perc=numeric())


cl1=6
cl2=5
#testa sobreposição de clusteres
for(cl1 in nxt[1]:top[1]){
  A<-pBreaks[[1]]$ini[cl1]
  B<-pBreaks[[1]]$fim[cl1]
  for(cl2 in nxt[2]:top[2]){
    C<-pBreaks[[2]]$ini[cl2]
    D<-pBreaks[[2]]$fim[cl2]
    #superpostos
    if(A<=D & B>=C){
      #AB contém em CD
      if(C>=A & D<=B){
        reg1<-c(cl1,cl2,2,1)
        reg2<-c(cl2,cl1,1,1)
        obs<-paste0("O cluster ", cl1," do intervalo 1 contém o cluster ", cl2," do intervalo 2.\n")
        }
      #AB contido em CD
      if(A>=C & B<=D){
        reg1<-c(cl1,cl2,1,1)
        reg2<-c(cl2,cl1,2,1)
        obs<-paste0("O cluster ", cl1," do intervalo 1 está contido no cluster ", cl2," do intervalo 2.\n")
      }
      if(A>=C & !B<=D){
        #superposição de Head em AB
        #tamanho<-ifelse((B-A>=D-C),D-C,B-A)
        reg1<-c(cl1,cl2,3,((D-A)/(B-A)))
        reg2<-c(cl2,cl1,4,((D-A)/(D-C)))
        obs<-paste0("O cluster ", cl1," do intervalo 1 tem ",round(((D-A)/(B-A)),digits = 2),
                    " de seu início dentro do cluster ", cl2," do intervalo 2.\n")
        }
      if(!A>=C & B<=D){
        #superposição de Tail em AB
        #tamanho<-ifelse((B-A>=D-C),D-C,B-A)
        reg1<-c(cl1,cl2,4,((B-C)/(B-A)))
        reg2<-c(cl2,cl1,3,((B-C)/(D-C)))
        obs<-paste0("O cluster ", cl1," do intervalo 1 tem ",round(((B-C)/(B-A)),digits = 2),
                    " de sua cauda dentro do cluster ", cl2," do intervalo 2.\n")
      }
    #Não superpostos
    }else{
      reg1<-c(cl1,cl2,0,0)
      reg2<-c(cl2,cl1,0,0)
    }
    if(reg1[3]!=0){
      int1to2<-rbind(int1to2,reg1)
      int2to1<-rbind(int2to1,reg2)
      cat(file = "/home/clovis/Dropbox/Chumbo/clusters.txt",obs,append = T)
    }
    correspond<-rbind(correspond,reg1)
  }
}
colnames(int1to2)<-c("int1","int2","rel","perc")
colnames(int2to1)<-c("int2","int1","rel","perc")
colnames(correspond)<-c("int1","int2","rel","perc")
duplic<-int1to2$int1[duplicated(int1to2$int1)]
obs<-paste0("\nClusters divididos :",toString(duplic),"\n")
cat(file = "/home/clovis/Dropbox/Chumbo/clusters.txt",obs,append = T)
int1to2<-int1to2[!(int1to2$int1%in%duplic & int1to2$perc<0.25),]

#correspond<-
#############################
#Isso está no lugar errado
#ver depois
itv1<-  read.table(file = "/home/clovis/Dropbox/Chumbo/terms/topTermsInterval1.csv", 
                    sep="\t",
                   header = T)
itv2<-  read.table(file = "/home/clovis/Dropbox/Chumbo/terms/topTermsInterval2.csv", 
                   sep="\t",
                   header = T)

GOs<-unique(rbind(itv1[,1:2],itv2[,1:2]))

itv1<-merge(itv1,int1to2[1:2],by.x="ClusterNumber", by.y = "int1")

ambos<-merge(GOs,itv1[,c(2,5)],by="GO.ID",all=T)
colnames(ambos)<-c("GO.ID", "Term", "Intv1")
ambos<-merge(ambos,itv2[,c(1,4)],by="GO.ID",all=T)
colnames(ambos)<-c("GO.ID", "Term", "Intv1", "Intv2")
ambos<-unique(ambos)
#verifica se cluster 1 é circular
if(crc[1]==1){
  ambos$Intv1[ambos$Intv1==top[1]]<-1
}
if(crc[2]==1){
  ambos$Intv2[ambos$Intv1==top[1]]<-1
}
ambos$cluster[!is.na(ambos$Intv1)]<-ambos$Intv1[!is.na(ambos$Intv1)]
ambos$cluster[!is.na(ambos$Intv2)]<-ambos$Intv2[!is.na(ambos$Intv2)]
ambos$Intv1[!is.na(ambos$Intv1)]<-"X"
ambos$Intv2[!is.na(ambos$Intv2)]<-"X"
ambos$Intv1[is.na(ambos$Intv1)]<-" "
ambos$Intv2[is.na(ambos$Intv2)]<-" "
ambos<-ambos[order(-ambos$cluster,ambos$Intv1,ambos$Intv2,decreasing = T),]
ambos<-ambos[,c(1,2,5,3,4)]
write.csv(ambos,file = "/home/clovis/Dropbox/Chumbo/terms/resumo.csv",
          row.names = F)

#######################################################3
#############################
# a correspondencia eu fiz na mão
#############################
#sem boundery
# clusteres<-data.frame(T1=c(1,2,3,4,NA,5,6,7,8,9,10,11,NA,NA),
#                       T2= c(1,1,1,2,3 ,4,5,6,7,8,8 ,9 ,10,11))
# save(clusteres,file = "./Data/clusteres.RData")
clusteres<-data.frame(T1=c(1,2,NA,3,NA,NA,NA,NA,NA,4,5,NA,NA,6,7,7,NA,8,9,10,11,12,13,NA,14,NA,NA,NA,15),
                      T2= c(1,2,3,4,5,6,7,8,9,10,NA,11,12,13,14,15,16,17,18,19,20,20,20,21,22,23,24,25,26))
save(clusteres,file = "./Data/clusteres.RData")
