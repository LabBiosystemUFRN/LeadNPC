rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadTest/"
setwd(baseDir)
source(file = "bin/00base.R")



# Ler tabela de counts ----------------------------------------------------


figuras="figuras"
graficos="graficos"

load(file = "./Data/counts.RData")

#create biomart list of association
# library('biomaRt')
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# #searchFilters(mart = mart, pattern = "ensembl.*id")
# genes <- rownames(logCPM)
# 
# G_list <- getBM(filters= "ensembl_gene_id", 
#                 attributes= c("ensembl_gene_id","hgnc_symbol"),
#                 mart= mart,
#                 values = genes)
#save(G_list,file = "./Data/Glist.RData")

load(file = "./Data/Glist.RData")
load(file = "./Data/counts.RData")
logCPM<-as.data.frame(logCPM)
logCPM$ensg<-rownames(logCPM)

logCPM<-merge(logCPM,G_list,by.x="ensg",by.y="ensembl_gene_id")


#duplicate day 0 to treated group
day0<-pheno_data[pheno_data$Day==0,]
day0$Group<-"Lead30"
pheno_data<-rbind(pheno_data,day0)

# ,"VIM" Radial glial cell markers
#"MAPT", Não achei
#"CXCR4", tronco para neural

# ncp<-c("PAX6","MSI1","NES","NOTCH1","SOX1","SOX2","FUT4")
# neuron<-c("TH","NEUROD6","DCX", "RBFOX3","BCL11B","SOX11","SOX4")
#dopaminergic markers TH, ,"LMX1B"
ncp<-c("MSI1","NES","NOTCH1","SOX1")
neuron<-c("TH","NEUROD6","DCX", "RBFOX3","BCL11B","SOX11","SOX4","LMX1B")

#to create color palette
# library(wesanderson)
# color=c(wes_palette("Rushmore1")[c(1,3,5)],
#         #wes_palette("Royal2"),
#         wes_palette("Royal1")[c(1)],
#         "#0000FF","#00FF00","#FF8C00",
#         wes_palette("BottleRocket2")[1:4])
#colourpicker::colourPicker(11)
# ggplot(colors, aes(factor(gene), fill=color)) +
#   geom_bar() +
#   scale_fill_manual(values = colors$color)

color<-c("#E1BD6D","#228B22","#F0800F","#899DA4","#0000FF","#00FF00","#8A2BE2","#FAD510","#CB2314","#273046","#354823","#00CDCD")
colors<-data.frame(gene=c(ncp,neuron),
                   color=color,
                   stringsAsFactors = F)
xcolors <- setNames(colors$color, colors$gene)

logCPM$type[logCPM$hgnc_symbol%in%ncp]<-1
logCPM$type[logCPM$hgnc_symbol%in%neuron]<-2
ctlExp<-as.data.frame(logCPM[logCPM$type%in%c(1,2),c(2:56)])

vgroup=c("Control","Lead30")
# vbase<-list(ctlExp[,c("SRR3944315","hgnc_symbol")],
#         ctlExp[,c("SRR3944352","hgnc_symbol")])
vbase<-list(ctlExp[,c("SRR3944315","hgnc_symbol")],
            ctlExp[,c("SRR3944315","hgnc_symbol")])


library(tidyr)

treat=1
g<-list()
gIndex=1
#loop ----
for(treat in c(1,2)){
  group<-vgroup[treat]
  base<-vbase[[treat]]
  # PAX6, SOX1
  # NEUROD6,RBFOX3, DCX, 
  ctlExp2<-gather(ctlExp,key,value,-54:-55)
  ctlExp2$value<-as.numeric(ctlExp2$value)
  ctlExp2<-ctlExp2[ctlExp2$key%in%pheno_data$Run[pheno_data$Group ==group],]
  
  control<-pheno_data[pheno_data$Group == group,]
  control<-control[order(control$Day),]
  
  colnames(base)<-c("base","hgnc_symbol")
  
  ctlExp2<- merge(ctlExp2,control[,c("Run","Day")],
                  by.x = "key",
                  by.y= "Run")
  ctlExp2<-merge(ctlExp2,base, by="hgnc_symbol")
  ctlExp2$value<-(ctlExp2$value/ctlExp2$base)
  ctlExp2<-ctlExp2[order(ctlExp2$Day),]
  
  ctlExp2<-merge(ctlExp2,
                 colors, 
                 by.x = "hgnc_symbol",
                 by.y = "gene")
  ctlExp2$maxX<-26
  
  library(ggplot2)
  library(ggrepel)
  
  miny<-min(ctlExp2$value[ctlExp2$type==1])
  maxy<-max(ctlExp2$value[ctlExp2$type==1])
  timeLabels<-data.frame(label=c("Time-interval 1","Time-interval 2"),
                         x=c((11-3)/2+3,(26-11)/2+11),
                         y=rep((maxy-miny)*0.7+miny,2))
  #plot 1----
  g[[gIndex]]<-ggplot()+theme_bw()+
    ggtitle(paste("NPC -",ifelse(treat==1,"Control","Treated")))+
    xlab("Days")+
    ylab("Relative expression")+
    geom_rect(aes_(xmin=3,
                  xmax=11.5,
                  ymin=miny,
                  ymax=maxy),
              fill="#fbf2cc60")+
    geom_rect(aes_(xmin=11.5,
                  xmax=26,
                  ymin=miny,
                  ymax=maxy),
              fill="#cce3f260")+
    geom_text(data = timeLabels, 
              aes(x,y, label=label),
              col=alpha(colour = "black",alpha = 0.25))+
    geom_line(data=ctlExp2[ctlExp2$type==1&
                             ctlExp2$hgnc_symbol!="PAX6",], 
              aes(x=Day, 
                  y=value,
                  color=hgnc_symbol))+
    scale_color_manual(values = xcolors)+
    guides(colour = "none", linetype="none")+
    geom_text_repel(data = ctlExp2[ctlExp2$type==1&
                                     ctlExp2$Day==26,],
                    aes(x=Day, 
                        y=value,
                        label=hgnc_symbol,
                        #hjust = "left",
                        color=hgnc_symbol),
                    segment.size = .1,
                    nudge_x = .1 ,
                    cex=3)
  
  gIndex<-gIndex+1
  
  miny<-min(ctlExp2$value[ctlExp2$type==2])
  maxy<-max(ctlExp2$value[ctlExp2$type==2])
  timeLabels<-data.frame(label=c("Time-interval 1","Time-interval 2"),
                         x=c((11-3)/2+3,(26-11)/2+11),
                         y=rep((maxy-miny)*0.7+miny,2))
  textLabels<-ctlExp2[ctlExp2$type==2&
                        ctlExp2$Day==26,
                      c("hgnc_symbol","value")]
  textLabels$x<-26
  
  #plot 2 ----
  g[[gIndex]]<-ggplot()+theme_bw()+
    ggtitle(paste("Neuron -",ifelse(treat==1,"Control","Treated")))+
    xlab("Days")+
    ylab("Relative expression")+
    geom_rect(aes_(xmin=3,
                  xmax=11.5,
                  ymin=miny,
                  ymax=maxy),
              fill="#fbf2cc60")+
    geom_rect(aes_(xmin=11.5,
                  xmax=26,
                  ymin=min(ctlExp2$value[ctlExp2$type==2]),
                  ymax=max(ctlExp2$value[ctlExp2$type==2])),
              fill="#cce3f260")+
    geom_text(data = timeLabels, 
              aes(x,y, label=label),
              col=alpha(colour = "black",alpha = 0.25))+
    geom_line(data=ctlExp2[ctlExp2$type==2,], 
              aes(x=Day, 
                  y=value,
                  color=hgnc_symbol))+
    scale_color_manual(values = xcolors)+
    guides(colour = "none", linetype="none")+
    geom_text_repel(data = ctlExp2[ctlExp2$type==2&
                                     ctlExp2$Day==26,],
                    aes(x=Day, 
                        y=value,
                        label=hgnc_symbol,
                        #hjust = "left",
                        color=hgnc_symbol),
                    segment.size = .1,
                    nudge_x = .1 ,
                    cex=3)
  
  gIndex<-gIndex+1
  
  library(dplyr)
  
  #sem grupo 0 ----
  ctlExp2$TI<-ifelse(ctlExp2$Day%in%(0:11),1,2)
  ctlTime<-ctlExp2[ctlExp2$Day>=0,]
  
  genes<-unique(ctlTime$hgnc_symbol)
  i="NOTCH1"
  result<-data.frame(gene=character(),
                     pval=numeric(),
                     type=numeric(),
                     stringsAsFactors = F)
  for (i in genes) {
    tmp1<-ctlTime$value[ctlTime$hgnc_symbol==i&
                          ctlTime$TI==1]
    tmp2<-ctlTime$value[ctlTime$hgnc_symbol==i&
                          ctlTime$TI==2]
    result[nrow(result)+1,"gene"]<-i
    result[nrow(result),"pval"]<-t.test(tmp1,tmp2)$p.value
    result[nrow(result),"type"]<-unique(ctlTime$type[ctlTime$hgnc_symbol==i])
  }
  result$pCor<-p.adjust(result$pval,method = "BH")
  result[result$pCor>0.01,]
  
  ctTime<-ctlTime%>%
    group_by(TI,type,hgnc_symbol)%>%
    summarise(value=mean(value))
  ctTime$type[ctTime$type==1]<-"NPC"
  ctTime$type[ctTime$type==2]<-"Neuron"
  
  ctTime<- merge(ctTime,result, by.x ="hgnc_symbol",
                 by.y = "gene")
  ctTime$line<-1
  ctTime$line[ctTime$pCor>0.01]<-2
  
  #plot 3 ----
  g[[gIndex]]<-ggplot()+theme_bw()+
    ggtitle(ifelse(treat==1,"Control","Treated"))+
    xlab("Time Interval")+
    ylab("Mean expression")+
    facet_grid(.~type.x)+
    geom_point(data=ctTime, 
               aes(x=as.factor(TI), 
                   y=value,
                   group=hgnc_symbol,
                   color=factor(hgnc_symbol)))+
    geom_line(data=ctTime, 
              aes(x=as.factor(TI), 
                  y=value,
                  group=hgnc_symbol,
                  linetype=as.factor(line),
                  color=factor(hgnc_symbol)))+
    geom_text_repel(data = ctTime[ctTime$TI==2,],
                    aes(x=as.factor(TI), 
                        y=value,
                        label=hgnc_symbol,
                        hjust = "left",
                        color=factor(hgnc_symbol)),
                    segment.size = .1,
                    nudge_x = .1  )+
    scale_y_log10()+
    guides(colour = "none", linetype="none")+
    scale_color_manual(values = xcolors)
    
  gIndex<-gIndex+1
  
  # com grupo 0 ----
  ctlExp2$TI<-0
  ctlExp2$TI[ctlExp2$Day%in%(3:11)]<-1
  ctlExp2$TI[ctlExp2$Day%in%(12:26)]<-2
  ctlTime<-ctlExp2[ctlExp2$Day>=0,]
  
  genes<-unique(ctlTime$hgnc_symbol)
  i="SOX11"
  result<-data.frame(gene=character(),
                     pval1=numeric(),
                     pval2=numeric(),
                     type=numeric(),
                     stringsAsFactors = F)
  for (i in genes) {
    tmp<-ctlTime[ctlTime$hgnc_symbol==i,
                         c("TI","value")]
    tmp$TI<-paste0("T",tmp$TI)
    anova(ctlTime)
    aov(value ~ TI, data = tmp)
    pairwise.t.test(tmp$value,tmp$TI,p.adjust.method = "BH")
    summary(TukeyHSD(aov(value ~ TI, data = tmp)))
    
    tmp0<-ctlTime$value[ctlTime$hgnc_symbol==i&
                          ctlTime$TI==0]
    tmp1<-ctlTime$value[ctlTime$hgnc_symbol==i&
                          ctlTime$TI==1]
    tmp2<-ctlTime$value[ctlTime$hgnc_symbol==i&
                          ctlTime$TI==2]
    result[nrow(result)+1,"gene"]<-i
    result[nrow(result),"pval1"]<-t.test(tmp0,tmp1)$p.value
    result[nrow(result),"pval2"]<-t.test(tmp1,tmp2)$p.value
    result[nrow(result),"type"]<-unique(ctlTime$type[ctlTime$hgnc_symbol==i])
  }
  result$pCor1<-p.adjust(result$pval1,method = "BH")
  result$pCor2<-p.adjust(result$pval2,method = "BH")
  result[result$pCor>0.01,]
  
  ctTime<-ctlTime%>%
    group_by(TI,type,hgnc_symbol)%>%
    summarise(value=mean(value))
  ctTime$type[ctTime$type==1]<-"NPC"
  ctTime$type[ctTime$type==2]<-"Neuron"
  
  ctTime<- merge(ctTime,result, by.x ="hgnc_symbol",
                 by.y = "gene")
  ctTime$line1<-1
  ctTime$line1[ctTime$pCor1>0.01]<-2
  ctTime$line2<-1
  ctTime$line2[ctTime$pCor2>0.01]<-2
  
  #plot 4 ----
  g[[gIndex]]<-ggplot()+theme_bw()+
    ggtitle(ifelse(treat==1,"Control","Treated"))+
    xlab("Time Interval")+
    ylab("Mean expression")+
    facet_grid(.~type.x)+
    geom_point(data=ctTime, 
               aes(x=as.factor(TI), 
                   y=value,
                   group=hgnc_symbol,
                   color=factor(hgnc_symbol)))+
    geom_line(data=ctTime[ctTime$TI%in%c(0:1),], 
              aes(x=as.factor(TI), 
                  y=value,
                  group=hgnc_symbol,
                  linetype=as.factor(line1),
                  color=factor(hgnc_symbol)))+
    geom_line(data=ctTime[ctTime$TI%in%c(1:2),], 
              aes(x=as.factor(TI), 
                  y=value,
                  group=hgnc_symbol,
                  linetype=as.factor(line2),
                  color=factor(hgnc_symbol)))+
    geom_text_repel(data = ctTime[ctTime$TI==2,],
                    aes(x=as.factor(TI), 
                        y=value,
                        label=hgnc_symbol,
                        hjust = "left",
                        color=factor(hgnc_symbol)),
                    segment.size = .1,
                    nudge_x = .1  )+
    scale_y_log10()+
    guides(colour = "none", linetype="none")+
    scale_color_manual(values = xcolors)
  
  
  gIndex<-gIndex+1
  
}

library(gridExtra)
library(grid)

lay <- rbind(c(1,2),
             c(3,4))

grid<-grid.arrange(g[[1]],g[[2]],g[[5]], g[[6]],layout_matrix = lay)

ggsave(filename = "/home/clovis/Dropbox/Chumbo/Marquers/markSup.pdf", 
         plot = grid, 
         device = "pdf", 
        
         scale = 2.5, 
         #width = 6.72, height = 2.98, units = "in",
         dpi = 300)
ggsave(filename = "/home/clovis/Dropbox/Chumbo/Marquers/markControl2T.pdf", 
       plot = g[[3]], 
       device = "pdf", 
       
       scale = 1.5, 
       #width = 6.72, height = 2.98, units = "in",
       dpi = 300)

ggsave(filename = "/home/clovis/Dropbox/Chumbo/Marquers/markControl3T.pdf", 
       plot = g[[4]], 
       device = "pdf", 
       
       scale = 1.5, 
       #width = 6.72, height = 2.98, units = "in",
       dpi = 300)

ggsave(filename = "/home/clovis/Dropbox/Chumbo/Marquers/markTreat2T.pdf", 
       plot = g[[7]], 
       device = "pdf", 
       
       scale = 1.5, 
       #width = 6.72, height = 2.98, units = "in",
       dpi = 300)

ggsave(filename = "/home/clovis/Dropbox/Chumbo/Marquers/markTreat3T.pdf", 
       plot = g[[8]], 
       device = "pdf", 
       
       scale = 1.5, 
       #width = 6.72, height = 2.98, units = "in",
       dpi = 300)
