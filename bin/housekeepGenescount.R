rm(list = ls())
library("edgeR")
library(biomaRt)
library(ggplot2)

setwd("/home/clovis/Dropbox/Chumbo/")


load(file = "./Data/counts.RData")
#load(file = "./Data/associationHs700.RData")

#genes house keeping beta-actin, GAPDH, G6PD
#beta-Actin and GAPDH housekeeping gene expression in asthmatic airways is variable and not suitable for normalising mRNA levels.
hkGenes<-c("ENSG00000075624","ENSG00000111640","ENSG00000160211")
exp<-as.data.frame(logCPM[rownames(logCPM)%in%hkGenes,])
controle<-as.data.frame(t(exp[,colnames(exp)%in%pheno_data$Run[pheno_data$Group=="Control"]]))
controle$Run <- rownames(controle)
controle<-merge(controle,pheno_data[,c(6,28)], by="Run")
caso<-as.data.frame(t(exp[,colnames(exp)%in%pheno_data$Run[pheno_data$Group=="Lead30"]]))
caso$Run <- rownames(caso)
caso<-merge(caso,pheno_data[,c(6,28)], by="Run")
lm1<-(rbind(controle[,c(2,5)], caso[,c(2,5)]))
lm2<-(rbind(controle[,c(3,5)], caso[,c(3,5)]))
lm3<-(rbind(controle[,c(4,5)], caso[,c(4,5)]))

sc1 <- as.data.frame(spline(controle$Day, controle$ENSG00000075624, n=100*length(controle)))
sc2 <- as.data.frame(spline(controle$Day, controle$ENSG00000111640, n=100*length(controle)))
sc3 <- as.data.frame(spline(controle$Day, controle$ENSG00000160211, n=100*length(controle)))
sl1 <- as.data.frame(spline(caso$Day, caso$ENSG00000075624, n=100*length(caso)))
sl2 <- as.data.frame(spline(caso$Day, caso$ENSG00000111640, n=100*length(caso)))
sl3 <- as.data.frame(spline(caso$Day, caso$ENSG00000160211, n=100*length(caso)))


p<- ggplot()+
  # stat_smooth(data=controle,aes(x=Day,y=ENSG00000075624, color = "1"),
  #             method = lm,
  #             formula = y ~ poly(x, length(controle)),se = T)+
  geom_smooth(data=lm1,aes(Day,ENSG00000075624,linetype="2"),
              color = "black",
              method = lm,
              se=F,
              size =0.2)+
  geom_smooth(data=lm2,aes(Day,ENSG00000160211,linetype="2"),
              color = "black",
              method = lm,
              se=F,
              size =0.2)+
  geom_smooth(data=lm3,aes(Day,ENSG00000111640,linetype="2"),
              color = "black",
              method = lm,
              se=F,
              size =0.2)+
geom_line(data=sc1,aes(x,y, color = "1"))+
  geom_point(data=controle,aes(Day,ENSG00000075624, shape ="1", color = "1"))+
  #geom_line(data=controle,aes(Day,ENSG00000111640, color = "2"))+
  geom_line(data=sc2,aes(x,y, color = "2"))+
  geom_point(data=controle,aes(Day,ENSG00000111640, shape ="1", color = "2"))+
  #geom_line(data=controle,aes(Day,ENSG00000160211, color = "3"))+
  geom_line(data=sc3,aes(x,y, color = "3"))+
  geom_point(data=controle,aes(Day,ENSG00000160211, shape ="1", color = "3"))+
  #geom_line(data=caso,aes(Day,ENSG00000075624, color = "1"))+
  geom_line(data=sl1,aes(x,y, color = "1"))+
  geom_point(data=caso,aes(Day,ENSG00000075624, shape ="2", color = "1"))+
  #geom_line(data=caso,aes(Day,ENSG00000111640, color = "2"))+
  geom_line(data=sl2,aes(x,y, color = "2"))+
  geom_point(data=caso,aes(Day,ENSG00000111640, shape ="2", color = "2"))+
  #geom_line(data=caso,aes(Day,ENSG00000160211, color = "3"))+
  geom_line(data=sl3,aes(x,y, color = "3"))+
  geom_point(data=caso,aes(Day,ENSG00000160211, shape ="2", color = "3"))+
  ggtitle("Housekeeping gene's expression through time")+
  xlab("Days")+
  ylab("LogCPM")+
  theme_bw()+
  scale_color_manual(name="Gene",
                     values=c("1"="red","2"="blue","3"="green"),
                     labels=c("ACTB", "GAPDH", "G6PD"," "))+
  scale_linetype_manual(name=" ",
                        values=c("2"= 2),
                        labels=c("Regression"))+
  scale_shape_manual(name="Sample",
                     values=c("1"=1,"2"=2),
                     labels=c("Control", "Lead30"))+
  scale_y_continuous(limits = c(0,max(exp)),breaks = seq(0,12,2))+
  scale_x_continuous(limits = c(0,26),breaks = seq(0,26,5))

pdf(width = 11,height = 8,file = paste0("./figuras/HKGenes.pdf"))
  suppressMessages(graphics::plot(p))
dev.off()
