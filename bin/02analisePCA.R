rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadTest/"
setwd(baseDir)
source(file = "bin/00base.R")

library(transcriptogramer)
library(biomaRt)
library("FactoMineR")
library("factoextra")
library(ggplot2)
library(factoextra)

# Ler tabela de counts ----------------------------------------------------


figuras="figuras"
graficos="tmpGraf"
transcripto="allTranscriptogramers80"
boundary=TRUE

load(file = "./Data/counts.RData")
load(file = "./Data/associationHs700.RData")

t(apply(association[1:100,], 1, sort))

# var set ----
name = "Lead30"
set<-"All"
radius<-80
pval <- 0.001
transc<-list()

#PCA ----
# PCAData0 <- transcriptogramPreprocess(association = "GenesS700.txt", ordering = "ordering_GenesS700.txt",
#                                 radius = 50 ) 
PCAData0 <- transcriptogramPreprocess(association = association, ordering = Hs700,
                                      radius = radius )#all genes

PCAData0 <- transcriptogramStep1(object = PCAData0, expression = logCPM,
                                 dictionary = dic_transcriptogramer, nCores = T)

g30 <- pheno_data$Run[pheno_data$source_name == "Control_NPCs" ]
PCAData<- PCAData0@transcriptogramS1[,colnames(PCAData0@transcriptogramS1)%in%g30]


g30 <- pheno_data$Run[pheno_data$source_name == "Lead30_NPCs"]
PCAData<- cbind(PCAData,PCAData0@transcriptogramS1[,colnames(PCAData0@transcriptogramS1)%in%g30])


names<-c(paste0(pheno_data$Day[pheno_data$source_name == "Control_NPCs"]),
         paste0(pheno_data$Day[pheno_data$source_name == "Lead30_NPCs"]))
# Nomear as colunas da tabela de expressao de acordo com os grupos
colnames(PCAData)<-names

#PCAData$Protein<-NUL/home/clovis/Dropbox/Clovis/aL
#PCAData$Position<-NULL
colnames(PCAData)<-names

targets <- as.factor(c(rep("Ctl", 27),
                       rep("L30", 26)))

# Estabelecer o esquema de cores
color.code<-colorRampPalette(c('blue','red'))(2)
pch.code<-c(15,16)


tPCAData<-t(PCAData)

# calcular a PCA
pca <- prcomp(tPCAData,scale = T)

# Plotar

pcs<-data.frame(pca$x)

shapes<-as.character(pch.code[targets])

g<-ggplot(data = pcs[,1:2],
       aes(PC1,PC2,
           label=names,
           shape=shapes,
           color = color.code[targets]))+
  geom_point()+
  geom_text(check_overlap = TRUE, cex=3, nudge_x = -3, col=1)+
  theme_bw()+
  scale_shape_manual(name = "",
                     guide = "legend",
                     labels=c("Control", "Lead 30"),
                     values = c(15,16))+
  scale_color_manual(name = "",
                     guide = "legend",
                     labels=c("Control", "Lead 30"),
                     values = color.code)+
  theme(legend.position = c(.1,0.2),
        legend.background = element_rect(fill = alpha("white",0)))

pdf(file = paste0("./",figuras,"/PC1xPC2.pdf"),height = 8, width = 11)
g
dev.off()


#clusters caso

tPCAData<-t(PCAData[,28:53])

res.pca<-PCA(tPCAData,ncp = 30, graph = F)

res.hcpc<-HCPC(res.pca,graph=T,nb.clust = -1,description = F)
#inverte
res.hcpc$call$X[,1:2]<- -(res.hcpc$call$X[,1:2])

fviz_dend(res.hcpc,
          cex = 0.7,
          palette = "jco",
          rect = T,
          rect_fill = T,
          rect_border = "jco",
          labels_track_height = 1,type = "phylogenic", repel = T, ggtheme = theme_bw())
dev.copy(pdf,paste0("./",figuras,"/Dendogram.pdf"))
dev.off()


g<-fviz_cluster(res.hcpc,
             cex = 0.7,
             axes = c(1,2),
             palette = "jco",
             rect = T,
             rect_fill = T,
             rect_border = "jco",
             labels_track_height = 1,
             type = "phylogenic", 
             repel = T, 
             ggtheme = theme_bw())
g<-g+  theme(legend.position = c(.9,0.2),
             legend.background = element_rect(fill = alpha("white",0)))+
  ggtitle("") +
  xlab("PC1") + ylab("PC2")


pdf(file = paste0("./",figuras,"/clusters.pdf"),height = 8, width = 11)
g
dev.off()



g<-fviz_screeplot(res.pca,ncp=17, 
                  geom = "bar",
                  ggtheme = theme_bw(),
                  choice = c("variance", "eigenvalue"))
sd<-as.data.frame(res.pca$eig)
#sd<-as.data.frame(t(sd$importance))
sd$pc<-c(1:nrow(sd))
sd<-rbind(c(0,0,0,0),sd)
colnames(sd)<-c("sd","prop","cum","pc")


g<-g+geom_line(data = sd[1:18,],
               aes(x = pc,y=cum,color="c"),
               lty=2)+
  ylim(0,100)+
  scale_color_manual(name="",
                     labels=c("Cumulative"),
                     values = c("c"="coral"))+
  theme(legend.position = c(.9,.5),
        legend.background = element_rect(fill = alpha("white",0)))+
  labs(title="")
pdf(file = paste0("./",figuras,"/variance.pdf"),height = 8, width = 11)
g
dev.off()



# CONTROL - LEAD3 ---------------------------------------------------------
# Criar objeto do transcriptogramer ---------------------------------------


# Passo 1 do transcriptogramer --------------------------------------------

i=2

for(i in seq(1,2)){
  if(i == 1){
    group<-c(3,4,5,6,7,8,9,10,11)
  }else
    if(i == 2){
      group<-c(12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
      }
  
  case_control30 <- pheno_data$Run[pheno_data$source_name == "Control_NPCs"&
                                     pheno_data$Day%in%group]
  
  
  
  # case_control30 <-append(case_control30, pheno_data$Run[pheno_data$source_name == paste0(name,"_NPCs")&
  #                                                          pheno_data$Day%in%c(i:(i+top-1))])
  case_control30 <- append(case_control30, pheno_data$Run[pheno_data$source_name == paste0(name,"_NPCs")&
                                                            pheno_data$Day%in%group])
  expression<- as.data.frame(logCPM[, colnames(logCPM)%in%case_control30])
  
  
  
  # t1 <- transcriptogramPreprocess(association = association3000S900, ordering = "ordering_genes3000S900.txt",
  #                                 radius = 25 )
  # t1 <- transcriptogramPreprocess(association = "/home/clovis/Doutorado/Artigos/Chumbo/GenesS700.txt", 
  #                                 ordering = "/home/clovis/Doutorado/Artigos/Chumbo/ordering_GenesS700.txt",
  #                                 radius = 50 )
  t1 <- transcriptogramPreprocess(association = association, 
                                  ordering = Hs700,
                                  radius = radius )
  
  t1 <- transcriptogramStep1(object = t1, expression = expression,
                             dictionary = dic_transcriptogramer, nCores = T)
  
  t1 <- transcriptogramStep2(object = t1, nCores = T)
  
  write.table(t1@transcriptogramS2, file = paste0("./samples/",set,"W",radius,"PCAGr",i,name,".csv"), sep="\t")
  
  
  levels <- c(rep(TRUE,length(group)),rep(FALSE,length(group)))
  
  levels <- c(rep(TRUE,length(group)),rep(FALSE,length(group)))
  
  tmp2<-rbind(levels,case_control30)
  write.table(tmp2, file = paste0("./levels/",set,"W",radius,"LevelsPCAGr",i,name,".csv"), sep="\t")
  
  
  possibleError<-tryCatch(t1 <- differentiallyExpressed(object = t1, levels = levels, pValue = pval, species = "Homo sapiens",
                                                        boundaryConditions = boundary,
                                                        title = paste0("Differential Expression - " ,name,"  - Group ",i," - Radius ",radius," - p-value ",pval)),
                          error=function(e) e)
  # possibleError<-tryCatch(t1 <- differentiallyExpressed(object = t1, levels = levels, pValue = 0.01,
  #                                                       title = paste0("Differential Expression - " ,name,"  - Day ",i,"-",(i+top-1) )),
  #                         error=function(e) e)
  
  if(inherits(possibleError, "error")) {
    cat(paste("Group ", i," - FAIL - Nothing differentially expressed") ,
        file=paste0("./levels/byGroupsPCA",set,"W",radius,"PCAGr",name,".log"),append=TRUE,sep="\n")
    next
  }
  dev.copy(pdf,paste0("./",graficos,"/",set,"W",radius,"PCAGr",i,"_",name,".pdf"))
  dev.off()
  
  
  # Enriquecimento ----------------------------------------------------------
  t1 <- clusterEnrichment(object = t1, species = "homo sapiens",
                              pValue = pval, nCores = T, onlyGenesInDE = F, algorithm = "parentchild")
  
  write.table(Terms(t1), file = paste0("./terms/",set,"W",radius,"PCAGr",i,name,".csv"), sep="\t")
  
  transc[[i]]<-t1
}

save(transc, file = paste0("./Data/allTranscriptogramers",radius))
