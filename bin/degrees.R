library(igraph)
nodesP<-subset(clusters12,Clust2==1,select = c("Protein") )
nodesP$Protein<-as.character(nodesP$Protein)

nodesF1<-subset(clusters12,Clust1==1,select = c("Protein") )
nodesF1$Protein<-as.character(nodesF1$Protein)

nodesF2<-subset(clusters12,Clust1==2,select = c("Protein") )
nodesF2$Protein<-as.character(nodesF2$Protein)

nodesF3<-subset(clusters12,Clust1==3,select = c("Protein") )
nodesF3$Protein<-as.character(nodesF3$Protein)

pai<-induced_subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesP$Protein])
F1<-induced_subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesF1$Protein])
F2<-induced_subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesF2$Protein])
F3<-induced_subgraph(graph=gMaster, v=V(gMaster)[V(gMaster)$name %in% nodesF3$Protein])


c0 <- degree_distribution(pai)
c1 <- degree_distribution(F1)
c2 <- degree_distribution(F2)
c3 <- degree_distribution(F3)
cl<-list()
comp<-max(length(c0),length(c1),length(c2),length(c3))
cl[[1]]<-data.frame(degree=seq(0,(length(c0)-1)),
                 p0=c0)
cl[[2]]<-data.frame(degree=seq(0,(length(c1)-1)),
                    p0=c1)
cl[[3]]<-data.frame(degree=seq(0,(length(c2)-1)),
                    p0=c2)
cl[[4]]<-data.frame(degree=seq(0,(length(c3)-1)),
                    p0=c3)

#resolve problemas de 0
#elimina posição com zero se a janela contiver menos que 5 zeros
for(turn in 1:4){
  clu2<-cl[[turn]]
  clu3<-clu2[c(1,nrow(clu2)),]
  i=300
  for(i in 2:(nrow(clu2)-1)){
    if(clu2[i,2]==0){
      if(sum(clu2[(i-1):(i+1),2])==0){
        clu3<-rbind(clu3,clu2[i,])
      }
    }else{
        clu3<-rbind(clu3,clu2[i,])
    }
  }
  cl[[turn]]<-clu3
}    
# smoothedLine<-stats::smooth.spline(clu2)
# df <- data.frame(x = smoothedLine$x, y = smoothedLine$y)
# 
# ggplot(df)+
#   geom_line(aes(x,y))+
#   scale_y_sqrt(limits = c(0,1))+
#   theme_bw()

colors<-c(1,2,3,4)
g<- ggplot()
for(i in 1:4){
  g<-g+geom_line(data=cl[[i]],aes(degree,p0), col=colors[i])
}
  g<-g+scale_y_sqrt(limits = c(0,1))+
  theme_bw()
g



require(graphics)
plot(dist ~ speed, data = cars, main = "data(cars)  &  smoothing splines")
cars.spl <- with(cars, smooth.spline(speed, dist))
cars.spl
## This example has duplicate points, so avoid cv = TRUE

lines(cars.spl, col = "blue")
ss10 <- smooth.spline(cars[,"speed"], cars[,"dist"], df = 10)
lines(ss10, lty = 2, col = "red")
legend(5,120,c(paste("default [C.V.] => df =",round(cars.spl$df,1)),
               "s( * , df = 10)"), col = c("blue","red"), lty = 1:2,
       bg = 'bisque')
