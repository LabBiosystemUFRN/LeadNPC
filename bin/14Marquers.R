rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadNPC/"
setwd(baseDir)
source(file = "bin/00base.R")


layout="col"
normalize = T



# Ler tabela de counts ----------------------------------------------------


figures="figures"
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
ncp<-c("MSI1","NES","SOX1","NOTCH1")
neuron<-c("TH","NEUROD6","DCX", "RBFOX3","GAD1","GAD2")

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

color<-c("#E1BD6D","#228B22","#F0800F","#899DA4","#0000FF","#00FF00","#8A2BE2","#FAD510","#CB2314","#00CDCD")
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
maxyMean<-0
maxyMean2<-0
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
  if(normalize){
    ctlExp2$value0<-(ctlExp2$value/ctlExp2$base)
  }else{  
    ctlExp2$value0<-ctlExp2$value
  }
  ctlExp2<-ctlExp2[order(ctlExp2$Day),]
  
  ctlExp2<-merge(ctlExp2,
                 colors, 
                 by.x = "hgnc_symbol",
                 by.y = "gene")
  ctlExp2$maxX<-26
  
  library(ggplot2)
  library(ggrepel)
  
  miny<-min(ctlExp2$value0[ctlExp2$type==1])
  maxy<-max(ctlExp2$value0[ctlExp2$type==1])
  if(normalize){
    timeLabels<-data.frame(label=c("Time-interval 1","Time-interval 2"),
                           x=c((11-3)/2+3,(26-11)/2+11),
                           y=rep((maxy-miny)*0.7+miny,2))
  }else{
    timeLabels<-data.frame(label=c("Time-interval 1","Time-interval 2"),
                           x=c((11-3)/2+3,(26-11)/2+11),
                           y=rep(1.25,2))
  }
  #plot 1----
  maxyMean2<-max(maxyMean2,max(ctlExp2$value0))
  
  g[[gIndex]]<-ggplot()+theme_bw()+
    ggtitle(paste("NPC -",ifelse(treat==1,"Control","Treated")))+
    xlab("Days")+
    ylab(ifelse(normalize,"Relative expression","Expression"))
  if(normalize){
    g[[gIndex]]<-g[[gIndex]]+
      geom_rect(aes_(xmin=3,
                     xmax=11.5,
                     ymin=miny,
                     ymax=maxy),
                fill="#fbf2cc60")+
      geom_rect(aes_(xmin=11.5,
                     xmax=26,
                     ymin=miny,
                     ymax=maxy),
                fill="#cce3f260")
  }else{
    g[[gIndex]]<-g[[gIndex]]+
      geom_rect(aes_(xmin=3,
                     xmax=11.5,
                     ymin=0,
                     ymax=maxyMean2),
                fill="#fbf2cc60")+
      geom_rect(aes_(xmin=11.5,
                     xmax=26,
                     ymin=0,
                     ymax=maxyMean2),
                fill="#cce3f260")+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1,
                                                        decimal.mark = '.'),
                         limits = c(0,maxyMean2))
  }
  g[[gIndex]]<-g[[gIndex]]+geom_text(data = timeLabels, 
            aes(x,y, label=label),
            col=alpha(colour = "black",alpha = 0.25))+
    geom_line(data=ctlExp2[ctlExp2$type==1&
                             ctlExp2$hgnc_symbol!="PAX6",], 
              aes(x=Day, 
                  y=value0,
                  color=hgnc_symbol))+
    scale_color_manual(values = xcolors)+
    guides(colour = "none", linetype="none")+
    geom_text_repel(data = ctlExp2[ctlExp2$type==1&
                                     ctlExp2$Day==26,],
                    aes(x=Day, 
                        y=value0,
                        label=hgnc_symbol,
                        #hjust = "left",
                        color=hgnc_symbol),
                    segment.size = .1,
                    nudge_x = .1 ,
                    cex=3)
  
  gIndex<-gIndex+1
  
  miny<-min(ctlExp2$value0[ctlExp2$type==2])
  maxy<-max(ctlExp2$value0[ctlExp2$type==2])
  
  if(normalize){
    timeLabels<-data.frame(label=c("Time-interval 1","Time-interval 2"),
                           x=c((11-3)/2+3,(26-11)/2+11),
                           y=rep((maxy-miny)*0.7+miny,2))
  }else{
    timeLabels<-data.frame(label=c("Time-interval 1","Time-interval 2"),
                           x=c((11-3)/2+3,(26-11)/2+11),
                           y=rep(1.25,2))
  }
  
  textLabels<-ctlExp2[ctlExp2$type==2&
                        ctlExp2$Day==26,
                      c("hgnc_symbol","value0")]
  textLabels$x<-26
  
  maxyMean2<-max(maxyMean2,max(ctlExp2$value0))
  
  #plot 2 ----
  g[[gIndex]]<-ggplot()+theme_bw()+
    ggtitle(paste("Neuron -",ifelse(treat==1,"Control","Treated")))+
    xlab("Days")+
    ylab(ifelse(normalize,"Relative expression","Expression"))+
    geom_rect(aes_(xmin=3,
                  xmax=11.5,
                  ymin=miny,
                  ymax=maxy),
              fill="#fbf2cc60")
  if(normalize){
    g[[gIndex]]<-g[[gIndex]]+geom_rect(aes_(xmin=11.5,
                                            xmax=26,
                                            ymin=miny,
                                            ymax=maxy),
                                       fill="#cce3f260")
  }else{
    g[[gIndex]]<-g[[gIndex]]+
      geom_rect(aes_(xmin=3,
                     xmax=11.5,
                     ymin=0,
                     ymax=maxyMean2),
                fill="#fbf2cc60")+
      geom_rect(aes_(xmin=11.5,
                     xmax=26,
                     ymin=0,
                     ymax=maxyMean2),
                fill="#cce3f260")+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1,
                                                        decimal.mark = '.'),
                         limits = c(0,maxyMean2))
  }
  g[[gIndex]]<-g[[gIndex]]+geom_text(data = timeLabels, 
                                     aes(x,y, label=label),
                                     col=alpha(colour = "black",alpha = 0.25))+
    geom_line(data=ctlExp2[ctlExp2$type==2,], 
              aes(x=Day, 
                  y=value0,
                  color=hgnc_symbol))+
    scale_color_manual(values = xcolors)+
    guides(colour = "none", linetype="none")+
    geom_text_repel(data = ctlExp2[ctlExp2$type==2&
                                     ctlExp2$Day==26,],
                    aes(x=Day, 
                        y=value0,
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
  
  maxyMean <-max(maxyMean ,max(ctTime$value))

  ctTime$type[ctTime$type==1]<-"NPC"
  ctTime$type[ctTime$type==2]<-"Neuron"
  
  ctTime<- merge(ctTime,result, by.x ="hgnc_symbol",
                 by.y = "gene")
  ctTime$line<-1
  ctTime$line[ctTime$pCor>0.01]<-2
  
  ctTime$type<-factor(ctTime$type.x,levels=c("NPC","Neuron"))
  
  #plot 3 ----
  g[[gIndex]]<-ggplot()+theme_bw()+
    ggtitle(ifelse(treat==1,"Control","Treated"))+
    xlab("Time Interval")+
    ylab("Mean expression")
  if(layout=="col"){
    g[[gIndex]]<-g[[gIndex]]+
      facet_grid(type~.)
  }else{
    g[[gIndex]]<-g[[gIndex]]+
      facet_grid(.~type)
  }
    g[[gIndex]]<-g[[gIndex]]+
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
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1,
                                                      decimal.mark = '.'),
                       limits = c(0,maxyMean))+
    guides(colour = "none", linetype="none")+
    scale_color_manual(values = xcolors)
    
  gIndex<-gIndex+1
  
  # com grupo 0 ----
  ctlExp2$TI<-0
  ctlExp2$TI[ctlExp2$Day%in%(3:11)]<-1
  ctlExp2$TI[ctlExp2$Day%in%(12:26)]<-2
  ctlTime<-ctlExp2[ctlExp2$Day>=0,]
  
  genes<-unique(ctlTime$hgnc_symbol)
  i="GAD2"
  result<-data.frame(gene=character(),
                     pval0to1=numeric(),
                     pval0to2=numeric(),
                     pval1to2=numeric(),
                     type=numeric(),
                     stringsAsFactors = F)
  for (i in genes) {
    tmp<-ctlTime[ctlTime$hgnc_symbol==i,
                         c("TI","value")]
    tmp$TI<-paste0("T",tmp$TI)
    #anova(ctlTime)
    #aov(value ~ TI, data = tmp)
    pairTest<-pairwise.t.test(tmp$value,tmp$TI,p.adjust.method = "BH")
    # summary(TukeyHSD(aov(value ~ TI, data = tmp)))
    # 
    # tmp0<-ctlTime$value[ctlTime$hgnc_symbol==i&
    #                       ctlTime$TI==0]
    # tmp1<-ctlTime$value[ctlTime$hgnc_symbol==i&
    #                       ctlTime$TI==1]
    # tmp2<-ctlTime$value[ctlTime$hgnc_symbol==i&
    #                       ctlTime$TI==2]
    # result[nrow(result),"pval1"]<-t.test(tmp0,tmp1)$p.value
    # result[nrow(result),"pval2"]<-t.test(tmp1,tmp2)$p.value
    # result[nrow(result),"pval3"]<-t.test(tmp0,tmp2)$p.value
    result[nrow(result)+1,"gene"]<-i
    result[nrow(result),"pval0to1"]<-pairTest$p.value[1,1]
    result[nrow(result),"pval0to2"]<-pairTest$p.value[2,1]
    result[nrow(result),"pval1to2"]<-pairTest$p.value[2,2]
    result[nrow(result),"type"]<-unique(ctlTime$type[ctlTime$hgnc_symbol==i])
  }
  # result$pCor1<-p.adjust(result$pval1,method = "BH")
  # result$pCor2<-p.adjust(result$pval2,method = "BH")
  # result[result$pCor>0.01,]
  # library(reshape2)
  # melt(result)
  
  ctTime<-ctlTime%>%
    group_by(TI,type,hgnc_symbol)%>%
    summarise(value=mean(value))
  
  maxyMean <-max(maxyMean ,max(ctTime$value))
  
  
  ctTime$type[ctTime$type==1]<-"NPC"
  ctTime$type[ctTime$type==2]<-"Neuron"
  
  ctTime<- merge(ctTime,result, by.x ="hgnc_symbol",
                 by.y = "gene")
  ctTime$line0to1<-1
  ctTime$line0to1[ctTime$pval0to1>0.01]<-2
  ctTime$line1to2<-1
  ctTime$line1to2[ctTime$pval1to2>0.01]<-2
  ctTime$line0to2<-1
  ctTime$line0to2[ctTime$pval0to2>0.01]<-0
  #plot 4 ----
  ctTime$type<-factor(ctTime$type.x,levels=c("NPC","Neuron"))
  
  g[[gIndex]]<-ggplot()+theme_bw()+
    ggtitle(ifelse(treat==1,"Control","Treated"))+
    xlab("Time Interval")+
    ylab("Mean expression")
  if(layout=="col"){
    g[[gIndex]]<-g[[gIndex]]+
      facet_grid(type~.)
  }else{
    g[[gIndex]]<-g[[gIndex]]+
      facet_grid(.~type)
  }
  g[[gIndex]]<-g[[gIndex]]+
    geom_point(data=ctTime, 
               aes(x=as.factor(TI), 
                   y=value,
                   group=hgnc_symbol,
                   color=factor(hgnc_symbol)),
               pch=20,cex = 1)+
    geom_line(data=ctTime[ctTime$TI%in%c(0:1),], 
              aes(x=as.factor(TI), 
                  y=value,
                  group=hgnc_symbol,
                  linetype=as.factor(line0to1),
                  color=factor(hgnc_symbol)))+
    geom_line(data=ctTime[ctTime$TI%in%c(1:2),], 
              aes(x=as.factor(TI), 
                  y=value,
                  group=hgnc_symbol,
                  linetype=as.factor(line1to2),
                  color=factor(hgnc_symbol)))+
    geom_point(data=ctTime[ctTime$TI%in%c(0,2)&
                             ctTime$line0to2==1,], 
              aes(x=as.factor(TI), 
                  y=value,
                  group=hgnc_symbol,
                  color=factor(hgnc_symbol)),
              pch = 15,cex=3)+
    geom_text_repel(data = ctTime[ctTime$TI==2,],
                    aes(x=as.factor(TI), 
                        y=value,
                        label=hgnc_symbol,
                        hjust = "left",
                        color=factor(hgnc_symbol)),
                    segment.size = .1,
                    nudge_x = .1,cex=6 )+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1,
                                                      decimal.mark = '.'),
                       limits = c(0,maxyMean))+
    guides(colour = "none", linetype="none")+
    scale_color_manual(values = xcolors)
  
  gIndex<-gIndex+1
  ctlExp2$treat<-group
  if(!exists("allTreat")){
    allTreat<-ctlExp2
  }else{
    allTreat<-rbind(allTreat,ctlExp2)
  }
}

#painel ----
#painel 1----

genes1<-unique(allTreat$hgnc_symbol[allTreat$type==1])
painel1<-list()
indice<-1
for(indice in 1:length(genes1)){
  maxy1<-max(allTreat$value0[allTreat$hgnc_symbol==genes1[indice]])
  kst<-allTreat[allTreat$hgnc_symbol==genes1[indice]&
                  allTreat$Day>0,]
  d1<-kst[kst$treat=="Control"&
            kst$Day%in%c(3:11),]
  d1<-d1$value0[order(d1$Day)]
  d2<-kst[kst$treat=="Lead30"&
            kst$Day%in%c(3:11),]
  d2<-d2$value0[order(d2$Day)]
  ks1<-ks.test(d1,d2, alternative = "t")

  d1<-kst[kst$treat=="Control"&
            kst$Day%in%c(12:26),]
  d1<-d1$value0[order(d1$Day)]
  d2<-kst[kst$treat=="Lead30"&
            kst$Day%in%c(12:26),]
  d2<-d2$value0[order(d2$Day)]
  ks2<-ks.test(d1,d2, alternative = "t")
  
  timeLabels<-data.frame(label=c("Time-interval 1","Time-interval 2"),
                         x=c((11-3)/2+3,(26-11)/2+11),
                         y=rep(maxy1/10,2))
  ks1<-p.adjust(c(ks1$p.value,ks2$p.value),
                method = "BH")
  pval1<-paste0("p-val: ",sprintf("%1.2E",ks1[1]))
  pval2<-paste0("p-val: ",sprintf("%1.2E",ks1[2]))

  painel1[[indice]]<-ggplot()+theme_bw()+
    ggtitle(genes1[indice])+
    xlab("Days")+
    ylab("Expression")+
    geom_rect(aes_(xmin=3,
                   xmax=11.5,
                   ymin=0,
                   ymax=maxy1),
              fill="#fbf2cc60")+
    geom_rect(aes_(xmin=11.5,
                   xmax=26,
                   ymin=0,
                   ymax=maxy1),
              fill="#cce3f260")+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1,
                                                      decimal.mark = '.'))+
    geom_text(data = timeLabels,
              aes(x,y, label=label),
              col=alpha(colour = "black",alpha = 0.25))+
    geom_text(aes_(x=7,y=(maxy1/4),
                   label=pval1),
              cex=3,
              col=ifelse(ks1[1]<=0.01,"blue","red"))+
    geom_text(aes_(x=19,y=(maxy1/4),
                   label=pval2),
              cex=3,
              col=ifelse(ks1[2]<=0.01,"blue","red"))+
    geom_line(data=allTreat[allTreat$hgnc_symbol==genes1[indice],], 
              aes(x=Day, 
                  y=value0,
                  color=hgnc_symbol,
                  linetype = factor(treat,levels=c("Lead30","Control"))))+
    # geom_line(data=allTreat[allTreat$treat == "Lead30"&
    #                           allTreat$hgnc_symbol==genes1[indice],], 
    #           aes(x=Day, 
    #               y=value0,
    #               color=hgnc_symbol),lty=1)+
    scale_color_manual(values = xcolors)+
    guides(colour = "none",linetype=guide_legend(" "))
  painel1[[1]]
}


library(gridExtra)
library(grid)

lay <- rbind(c(1,2),
             c(3,4))

grid<-grid.arrange(painel1[[1]],painel1[[2]],painel1[[3]], painel1[[4]],layout_matrix = lay)

ggsave(filename = paste0("/home/clovis/Dropbox/Chumbo/Marquers/FigS041a11",
                         ifelse(normalize,"Norm","Abs"),".pdf"), 
       plot = grid, 
       device = "pdf", 
       
       #scale = 2.5, 
       width = 14, height = 8.5, units = "in",
       dpi = 300)

#painel 2 ----

genes1<-unique(allTreat$hgnc_symbol[allTreat$type==2])
painel2<-list()
indice<-1
for(indice in 1:length(genes1)){
  maxy1<-max(allTreat$value0[allTreat$hgnc_symbol==genes1[indice]])
  kst<-allTreat[allTreat$hgnc_symbol==genes1[indice]&
                  allTreat$Day>0,]
  d1<-kst[kst$treat=="Control"&
            kst$Day%in%c(3:11),]
  d1<-d1$value0[order(d1$Day)]
  d2<-kst[kst$treat=="Lead30"&
            kst$Day%in%c(3:11),]
  d2<-d2$value0[order(d2$Day)]
  ks1<-ks.test(d1,d2, alternative = "t")
  
  d1<-kst[kst$treat=="Control"&
            kst$Day%in%c(12:26),]
  d1<-d1$value0[order(d1$Day)]
  d2<-kst[kst$treat=="Lead30"&
            kst$Day%in%c(12:26),]
  d2<-d2$value0[order(d2$Day)]
  ks2<-ks.test(d1,d2, alternative = "t")
  
  timeLabels<-data.frame(label=c("Time-interval 1","Time-interval 2"),
                         x=c((11-3)/2+3,(26-11)/2+11),
                         y=rep(maxy1/10,2))
  
  ks1<-p.adjust(c(ks1$p.value,ks2$p.value),
           method = "BH")
  pval1<-paste0("p-val: ",sprintf("%1.2E",ks1[1]))
  pval2<-paste0("p-val: ",sprintf("%1.2E",ks1[2]))
  
  painel2[[indice]]<-ggplot()+theme_bw()+
    ggtitle(genes1[indice])+
    xlab("Days")+
    ylab("Expression")+
    geom_rect(aes_(xmin=3,
                   xmax=11.5,
                   ymin=0,
                   ymax=maxy1),
              fill="#fbf2cc60")+
    geom_rect(aes_(xmin=11.5,
                   xmax=26,
                   ymin=0,
                   ymax=maxy1),
              fill="#cce3f260")+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1,
                                                      decimal.mark = '.'))+
    geom_text(data = timeLabels,
              aes(x,y, label=label),
              col=alpha(colour = "black",alpha = 0.25))+
    geom_text(aes_(x=7,y=(maxy1/4),
                   label=pval1),
              cex=3,
              col=ifelse(ks1[1]<=0.01,"blue","red"))+
    geom_text(aes_(x=19,y=(maxy1/4),
                   label=pval2),
              cex=3,
              col=ifelse(ks1[2]<=0.01,"blue","red"))+
    geom_line(data=allTreat[allTreat$hgnc_symbol==genes1[indice],], 
              aes(x=Day, 
                  y=value0,
                  color=hgnc_symbol,
                  linetype = factor(treat,levels=c("Lead30","Control"))))+
    # geom_line(data=allTreat[allTreat$treat == "Lead30"&
    #                           allTreat$hgnc_symbol==genes1[indice],], 
    #           aes(x=Day, 
    #               y=value0,
    #               color=hgnc_symbol),lty=1)+
    scale_color_manual(values = xcolors)+
    guides(colour = "none",linetype=guide_legend(" "))
}



lay <- rbind(c(1,2),
             c(3,4),
             c(5,6))

grid<-grid.arrange(painel2[[1]],painel2[[3]],
                   painel2[[2]], painel2[[4]],
                   painel2[[5]], painel2[[6]]
                   ,layout_matrix = lay)

ggsave(filename = paste0("/home/clovis/Dropbox/Chumbo/Marquers/FigS053a11",
                         ifelse(normalize,"Norm","Abs"),".pdf"), 
       plot = grid, 
       device = "pdf", 
       
       #scale = 2.5, 
       width = 14, height = 8.5, units = "in",
       dpi = 300)



library(gridExtra)
library(grid)
#save ----
#  g[[8]]
lay <- rbind(c(1,2),
             c(3,4))

grid<-grid.arrange(g[[1]],g[[2]],g[[5]], g[[6]],layout_matrix = lay)

ggsave(filename = paste0("./figures/markSup",
                         ifelse(normalize,"Norm","Abs"),".pdf"), 
         plot = grid, 
         device = "pdf", 
        
         #scale = 2.5, 
         width = 14, height = 8.5, units = "in",
         dpi = 300)

ggsave(filename = "/home/clovis/Dropbox/Chumbo/Marquers/markControl2T.pdf", 
       plot = g[[3]], 
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

       device = "pdf", 
       
       #scale = 1.8, 
       width = 11, height = 7, units = "in",
       dpi = 300)

