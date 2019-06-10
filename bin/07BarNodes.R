rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadNPC/"
setwd(baseDir)
source(file = "bin/00base.R")

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

load("./Data/allTranscriptogramers80")
load("./Data/tfs.RData")
load("./Data/colors.RData")

#figures="figures"

library(igraph)
library(RedeR)
library(transcriptogramer)
library(ggplot2)
library(purrr)
library(scales)

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
load(file = "./Data/clusteres.RData")
#load(file = "./Data/clusteresSB.RData")
c1<-clusteres$T1
c2<-clusteres$T2

load(file='./Data/clusters12.RData')
#load(file='./Data/clusters12SB.RData')

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

#maximo valor para o gráfico
maxval=max(table(cl12$Clust1[cl12$Clust1!=0]),
           table(cl12$Clust2[cl12$Clust2!=0]))
ref2=2
#setwd("/home/clovis/Dropbox/Chumbo/figures/Figura2/nodes/")
# clusters 1 p/ 1----
#somente quem é 1 para 1
for(ref2 in c(2,3,4,5,6,7,9,10,11)){
  ref1<-cltRef$c1[cltRef$c2 ==  ref2]
  ref1<-na.exclude(ref1)
  
  #Lista genes em cada tempo
  genes1<-nrow(cl12[cl12$Clust1%in%ref1,])
  genes2<-nrow(cl12[cl12$Clust2%in%ref2,])
  
  m<-data.frame(i=c(1,2),genes=(c(genes1/1000,genes2/1000)))
  
  theme_new <- theme_set(theme_bw(base_size = 50))
  
  theme_new <- theme_update(axis.title.x = element_text(size =30),
                            axis.title.y = element_text(size =30),
                            axis.text.x = element_text(size =50),
                            axis.text.y = element_text(size =50))
  
  g <- ggplot(m, aes(as.factor(i),genes))+
    #scale_x_continuous(limits = c(0.5,2.5))+
    scale_y_continuous(limits = c(0,maxval/1000))+
    xlab("Interval")+
    ylab("Number of nodes (x 1000)")+
    geom_bar(stat = "identity", col=1, fill =c(color[[2]][ref2],
                                               color[[2]][ref2]) ,
             width = 0.45)
  # geom_bar(stat = "identity", col=1, fill =c(color[[2]][ref2],
  #                                            lightColor(color[[2]][ref2])) ,
  #          width = 0.45)
  
  plot(g)
  
  svg(width = 11,height = 8.95,file = paste0("./figures/countGenesClust",ref2,".svg"))
  suppressMessages(graphics::plot(g))
  dev.off()
  
}
#clusteres 1 e 8
# cluster 1 ----
ref2<-1
ref11<-1
ref12<-2
ref13<-3

#Lista genes em cada tempo
genes11<-nrow(cl12[cl12$Clust1%in%ref11,])
genes12<-nrow(cl12[cl12$Clust1%in%ref12,])
genes13<-nrow(cl12[cl12$Clust1%in%ref13,])
genes2<-nrow(cl12[cl12$Clust2%in%ref2,])

m<-data.frame(i=c(1,1,1,2),
              genes=c(genes11/1000,
                      genes12/1000,
                      genes13/1000,
                      genes2/1000))

theme_new <- theme_set(theme_bw(base_size = 50))

theme_new <- theme_update(axis.title.x = element_text(size =30),
                          axis.title.y = element_text(size =30),
                          axis.text.x = element_text(size =50),
                          axis.text.y = element_text(size =50))

g <- ggplot()+
  scale_x_discrete(breaks = c(1,2))+
  scale_y_continuous(limits = c(0,maxval/1000))+
  xlab("Interval")+
  ylab("Number of nodes (x 1000)")+
  geom_bar(data=m[1:3,], aes(as.factor(i),genes),
           fill = c(color[[1]][ref11],color[[1]][ref12],color[[1]][ref13]),
           stat = "identity", col=1, width = 0.45)+
  geom_bar(data=m[4,], aes(as.factor(i),genes),
           fill = color[[2]][ref2],
           stat = "identity", col=1, width = 0.45)+
  # geom_bar(data=m[4,], aes(as.factor(i),genes),
  #          fill = lightColor(color[[2]][ref2]),
  #          stat = "identity", col=1, width = 0.45)+
  guides(fill=FALSE)

plot(g)

svg(width = 11,height = 8.95,file = paste0("./figures/countGenesClust",ref2,".svg"))
suppressMessages(graphics::plot(g))
dev.off()

# cluster 8 ----
ref2<-8
ref11<-9
ref12<-10

#Lista genes em cada tempo
genes11<-nrow(cl12[cl12$Clust1%in%ref11,])
genes12<-nrow(cl12[cl12$Clust1%in%ref12,])
genes2<-nrow(cl12[cl12$Clust2%in%ref2,])

m<-data.frame(i=c(1,1,2),
              genes=c(genes11/1000,
                      genes12/1000,
                      genes2/1000))

theme_new <- theme_set(theme_bw(base_size = 50))

theme_new <- theme_update(axis.title.x = element_text(size =30),
                          axis.title.y = element_text(size =30),
                          axis.text.x = element_text(size =50),
                          axis.text.y = element_text(size =50))

g <- ggplot()+
  scale_x_discrete(breaks = c(1,2))+
  scale_y_continuous(limits = c(0,maxval/1000))+
  xlab("Interval")+
  ylab("Number of nodes (x 1000)")+
  geom_bar(data=m[1:2,], aes(as.factor(i),genes),
           fill = c(color[[1]][ref11],color[[1]][ref12]),
           stat = "identity", col=1, width = 0.45)+
  geom_bar(data=m[3,], aes(as.factor(i),genes),
           fill = color[[2]][ref2],
           stat = "identity", col=1, width = 0.45)+
  # geom_bar(data=m[3,], aes(as.factor(i),genes),
  #          fill = lightColor(color[[2]][ref2]),
  #          stat = "identity", col=1, width = 0.45)+
  guides(fill=FALSE)

plot(g)

svg(width = 11,height = 8.95,file = paste0("./figures/countGenesClust",ref2,".svg"))
suppressMessages(graphics::plot(g))
dev.off()

#fim############################################################################  
