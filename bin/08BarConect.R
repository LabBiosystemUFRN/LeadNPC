rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadNPC/"
setwd(baseDir)
source(file = "bin/00base.R")

load("./Data/allTranscriptogramers80")
load("./Data/tfs.RData")
load("./Data/colors.RData")
#load("./Data/colorsSB.RData")
load("./Data/assocNoDup.RData")

conectDir="tmpGraf/"

library(igraph)
library(RedeR)
library(transcriptogramer)
library(ggplot2)
library(purrr)

#intersecção dos clusteres
 load(file = "./Data/clusteres.RData")
 c1<-clusteres$T1
 c2<-clusteres$T2

 load(file = "./Data/cl12.RData")
 load(file='./Data/clusters12.RData')

 #Grapho master ----

gMaster <- graph_from_data_frame(assocNoDup, directed = F)
V(gMaster)$nodeShape <- "CIRCLE"
E(gMaster)$edgeColor<-"#c9c9c9ff"

cltRef<-data.frame(c1=c1,c2=c2)

  # Intervalo 2 ----
  i=1
  netMean<-data.frame(time=numeric(), cl=numeric(), mean=numeric())
  for(i in 1:11){
    nodesP<-subset(clusters12,Clust2==i,select = c("Protein") )
    nodesP$Protein<-as.character(nodesP$Protein)
  
    F0<-induced_subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
    degree<-igraph::degree(F0, mode = "all", loops = F)
    netMean<-rbind(netMean,c(2,i,mean(degree)))
    c0 <- degree_distribution(F0)
    cl<-data.frame(degree=seq(0,(length(c0)-1)),
                      p0=c0)
    if(nrow(cl)<100){
      step<-1
    }else{
      step<-round(nrow(cl)/100,0)
    }

    clu3<-cl[c(1),]
    j=1
    for(j in seq(2,nrow(cl)-step,step)){
      clu3<-rbind(clu3,c(j,sum(cl[j:(j+step),2])))
    }
    clu3<-rbind(clu3,c(nrow(cl),sum(cl[j:nrow(cl),2])))
    clu3<-rbind(clu3,c(nrow(cl)+step,0))
    

    cl<-clu3
    cl$p0[cl$p0==0]<- 1e-15
    cl2<-data.frame(degree=log10(cl$degree),p0=log10(cl$p0))

    
    cl2$degree[cl2$degree==-Inf]<-0
    smoothedLine<-stats::smooth.spline(cl2)
    cl2 <- data.frame(degree = smoothedLine$x, p0 = smoothedLine$y)
    cl2$p0[cl2$p0 < -15]<- -15
    g<- ggplot()+
      geom_line(data=cl2,aes(degree,p0))+
      geom_hline(yintercept = -15,lty=3, col="red")+
      scale_y_continuous(breaks = c(-12.5,-10,-5,0))+
      xlab("Degrees of Connectivity (log)")+
      ylab("P(k) (log)")+
      annotate(geom="text", x=0, y=-15, label="-Inf",
               color="red",hjust = 3)+
      coord_cartesian(clip = 'off' ) +
      theme_bw()
    #    g
    svg(file = paste0(conectDir,"t2c",
                      i,".svg"),
        width = 11,height = 8)
    plot(g)
    dev.off()
}
  

  # Intervalo 1 ----
  i=1
  for(i in 1:11){
    nodesP<-subset(clusters12,Clust1==i,select = c("Protein") )
    nodesP$Protein<-as.character(nodesP$Protein)
    
    F0<-induced_subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
    c0 <- degree_distribution(F0)
    cl<-data.frame(degree=seq(0,(length(c0)-1)),
                   p0=c0)
    if(nrow(cl)<100){
      step<-1
    }else{
      step<-round(nrow(cl)/100,0)
    }
    #clu3<-cl[c(1,nrow(cl)),]
    clu3<-cl[c(1),]
    j=1
    for(j in seq(2,nrow(cl)-step,step)){
      clu3<-rbind(clu3,c(j,sum(cl[j:(j+step),2])))
    }
    clu3<-rbind(clu3,c(nrow(cl),sum(cl[j:nrow(cl),2])))
    clu3<-rbind(clu3,c(nrow(cl)+step,0))
    
    cl<-clu3
    cl1<-clu3
    cl$p0[cl$p0==0]<- 1e-15
    cl2<-data.frame(degree=log10(cl$degree),p0=log10(cl$p0))
    cl3<-data.frame(degree=log10(cl$degree),p0=cl$p0)

    cl2$degree[cl2$degree==-Inf]<-0
    cl3$degree[cl3$degree==-Inf]<-0
    
    smoothedLine<-stats::smooth.spline(cl1)
    cl1 <- data.frame(degree = smoothedLine$x, p0 = smoothedLine$y)

    smoothedLine<-stats::smooth.spline(cl2)
    cl2 <- data.frame(degree = smoothedLine$x, p0 = smoothedLine$y)
    cl2$p0[cl2$p0 < -15]<- -15

    smoothedLine<-stats::smooth.spline(cl3)
    cl3 <- data.frame(degree = smoothedLine$x, p0 = smoothedLine$y)

    g<- ggplot()+
      geom_line(data=cl3,aes(degree,p0))+
      xlab("Degrees of Connectivity (log)")+
      ylab("P(k)")+
      theme_bw()
    svg(file = paste0(conectDir,"t1c",
                      i,"logX.svg"),
        width = 11,height = 8)
    plot(g)
    dev.off()

    g<- ggplot()+
      geom_line(data=cl1,aes(degree,p0))+
      xlab("Degrees of Connectivity")+
      ylab("P(k)")+
      theme_bw()
    svg(file = paste0(conectDir,"t1c",
                      i,".svg"),
        width = 11,height = 8)
    plot(g)
    dev.off()
    
    g<- ggplot()+
      geom_line(data=cl2,aes(degree,p0))+
      geom_hline(yintercept = -15,lty=3, col="red")+
      scale_y_continuous(breaks = c(-12.5,-10,-5,0))+
      xlab("Degrees of Connectivity (log)")+
      ylab("P(k) (log)")+
      annotate(geom="text", x=0, y=-15, label="-Inf",
               color="red",hjust = 3)+
      coord_cartesian(clip = 'off' ) +
      theme_bw()
    svg(file = paste0(conectDir,"t1c",
                      i,"logXY.svg"),
        width = 11,height = 8)
    plot(g)
    dev.off()
    
    
  }
  
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
  
  
  
# conectiv média ----
  #rede pai
  cl=7
  netMean<-data.frame(time=numeric(), cl=numeric(), mean=numeric())
  time=1
  for(time in 1:2){
    clTmp<-na.exclude(clusters12[,c(1,time+2)])
    for(cl in 1:11){
      nodesP<-clTmp[clTmp[2]==cl,1]
      nodesP<-as.character(nodesP)
      # subAssoc<-association7[association7$p1%in%nodesP &
      #                         association7$p2%in%nodesP,]
      
      F0<-induced_subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP])
      degree<-igraph::degree(F0, mode = "out", loops = F)

      netMean[nrow(netMean)+1,] <- c(time,cl,mean(degree))
    }
  }
  
  rm(clTmp)
  conecDir="./figures/"
  ref2=3
  # clusters 1 p/ 1----
  #somente quem é 1 para 1
  for(ref2 in c(2,3,4,5,6,7,9,10,11)){
    if(nrow(na.exclude(cltRef[cltRef$c2 ==  ref2,]))==0){
      ref1=0
    }else{
      ref1<-cltRef$c1[cltRef$c2 ==  ref2]
      ref1<-na.exclude(ref1)
    }
    grData<-data.frame(cl=numeric(), media=numeric(),cor=character(),
                       stringsAsFactors = F)
    #parametros de cor tempo1
    cor<-substr(toupper(color[[1]][ref1]),1,7)
    if(ref1 == 0){
      grData[nrow(grData)+1,] <- data.frame(1,
                                            0,
                                            "#000000",
                                            stringsAsFactors = F)
    }else{
      grData[nrow(grData)+1,] <- data.frame(1,
                                 netMean$mean[netMean$time==1 &
                                                netMean$cl == ref1],
                                     cor,
                                 stringsAsFactors = F)
    }
    #parametros de cor tempo2
    cor<-toupper(color[[2]][ref2])
    # cor<-lightColor(toupper(color[[2]][ref2]))
    grData[nrow(grData)+1,] <- data.frame(2,
                                 netMean$mean[netMean$time==2 &
                                                netMean$cl == ref2],
                                 cor,
                                 stringsAsFactors = F)
    
    theme_new <- theme_set(theme_bw(base_size = 50))
    
    theme_new <- theme_update(axis.title.x = element_text(size =30),
                              axis.title.y = element_text(size =30),
                              axis.text.x = element_text(size =50),
                              axis.text.y = element_text(size =50))
    #grData$media2<-numeric(grData$media)
    #usar em caso de lightColor()
    #cor<-grData$cor
    cor<-c(cor,cor)
    g <- ggplot(grData, aes(x=as.factor(cl),y=media))+
      #scale_x_continuous(limits = c(0.5,2.5))+
      scale_y_continuous(limits=c(0,max(netMean$mean)))+
      xlab("Interval")+
      ylab("<k>")+
      geom_bar(aes(fill=cor),stat = "identity", col=1, width = 0.45)+
      scale_fill_manual(values=cor, guide=F)
    
    plot(g)
    
    svg(width = 11,height = 8.95,file = paste0(conecDir,"avgConec",ref2,".svg"))
    suppressMessages(graphics::plot(g))
    dev.off()
    
  }
  #clusteres 1 e 8
  # cluster 1 ----
  ref2<-1
  ref11<-1
  ref12<-2
  ref13<-3

  grData<-data.frame(cl=numeric(), media=numeric(),cor=character(),
                     stringsAsFactors = F)
  #parametros de cor tempo1
  grData[nrow(grData)+1,] <- data.frame(1,
                                        netMean$mean[netMean$time==1 &
                                                       netMean$cl == ref11],
                                        substr(toupper(color[[1]][ref11]),1,7),
                                          stringsAsFactors = F)
  grData[nrow(grData)+1,] <- data.frame(1,
                                        netMean$mean[netMean$time==1 &
                                                       netMean$cl == ref12],
                                        substr(toupper(color[[1]][ref12]),1,7),
                                        stringsAsFactors = F)
  grData[nrow(grData)+1,] <- data.frame(1,
                                        netMean$mean[netMean$time==1 &
                                                       netMean$cl == ref13],
                                        substr(toupper(color[[1]][ref13]),1,7),
                                        stringsAsFactors = F)
  #parametros de cor tempo2
  grData[nrow(grData)+1,] <- data.frame(2,
                                        netMean$mean[netMean$time==2 &
                                                       netMean$cl == ref2],
                                        # lightColor(toupper(color[[2]][ref2])),
                                        # stringsAsFactors = F)
                                        toupper(color[[2]][ref2]),
                                        stringsAsFactors = F)
  
  theme_new <- theme_set(theme_bw(base_size = 50))
  
  theme_new <- theme_update(axis.title.x = element_text(size =30),
                            axis.title.y = element_text(size =30),
                            axis.text.x = element_text(size =50),
                            axis.text.y = element_text(size =50))
  #grData$media2<-numeric(grData$media)
  cor<-grData$cor
  
  g <- ggplot()+
    scale_x_discrete(breaks = c(1,2))+
    scale_y_continuous(limits=c(0,max(netMean$mean)))+
    xlab("Interval")+
    ylab("<k>")+
    geom_bar(data=grData[1:3,], aes(x=as.factor(cl),y=media),
             fill = c(cor[1],cor[2],cor[3]),
             stat = "identity", col=1, width = 0.45)+
    geom_bar(data=grData[4,], aes(x=as.factor(cl),y=media),
             fill = c(cor[4]),
             stat = "identity", col=1, width = 0.45)
    #geom_bar(aes(fill=cor),stat = "identity", col=1, width = 0.45)+
    #scale_fill_manual(values=grData$cor, guide=F)
  
  plot(g)
  
  svg(width = 11,height = 8.95,file = paste0(conecDir,"avgConec1.svg"))
  suppressMessages(graphics::plot(g))
  dev.off()

  # cluster 8 ----
  ref2<-8
  ref11<-9
  ref12<-10
  
  grData<-data.frame(cl=numeric(), media=numeric(),cor=character(),
                     stringsAsFactors = F)
  #parametros de cor tempo1
  grData[nrow(grData)+1,] <- data.frame(1,
                                        netMean$mean[netMean$time==1 &
                                                       netMean$cl == ref11],
                                        substr(toupper(color[[1]][ref11]),1,7),
                                        stringsAsFactors = F)
  grData[nrow(grData)+1,] <- data.frame(1,
                                        netMean$mean[netMean$time==1 &
                                                       netMean$cl == ref12],
                                        substr(toupper(color[[1]][ref12]),1,7),
                                        stringsAsFactors = F)
  #parametros de cor tempo2
  grData[nrow(grData)+1,] <- data.frame(2,
                                        netMean$mean[netMean$time==2 &
                                                       netMean$cl == ref2],
                                        # lightColor(toupper(color[[2]][ref2])),
                                        # stringsAsFactors = F)
                                        toupper(color[[2]][ref2]),
                                        stringsAsFactors = F)
  
  theme_new <- theme_set(theme_bw(base_size = 50))
  
  theme_new <- theme_update(axis.title.x = element_text(size =30),
                            axis.title.y = element_text(size =30),
                            axis.text.x = element_text(size =50),
                            axis.text.y = element_text(size =50))
  #grData$media2<-numeric(grData$media)
  cor<-grData$cor
  
  g <- ggplot()+
    scale_x_discrete(breaks = c(1,2))+
    scale_y_continuous(limits=c(0,max(netMean$mean)))+
    xlab("Interval")+
    ylab("<k>")+
    geom_bar(data=grData[1:2,], aes(x=as.factor(cl),y=media),
             fill = c(cor[1],cor[2]),
             stat = "identity", col=1, width = 0.45)+
    geom_bar(data=grData[3,], aes(x=as.factor(cl),y=media),
             fill = c(cor[3]),
             stat = "identity", col=1, width = 0.45)
  #geom_bar(aes(fill=cor),stat = "identity", col=1, width = 0.45)+
  #scale_fill_manual(values=grData$cor, guide=F)
  
  plot(g)
  
  svg(width = 11,height = 8.95,file = paste0(conecDir,"avgConec8.svg"))
  suppressMessages(graphics::plot(g))
  dev.off()
  
  