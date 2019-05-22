
rm(list = ls())

#library(edgeR)
#library(refGenome)
library(transcriptogramer)
library(biomaRt)
library("FactoMineR")
library("factoextra")
library(ggplot2)
library(RedeR)
library(igraph)
library(dplyr)

library(data.table)
library(stringr)
library(limma)
library(edgeR)
library(readxl)
library(RColorBrewer)



# Ler tabela de counts ----------------------------------------------------
setwd("/home/clovis/Dropbox/Chumbo/Data/")
load("./allTranscriptogramers80")
load(file = "./counts.RData")
load(file = "./associationHs700.RData")
load(file = "./CPM.RData")


name = "Lead30"

#seleciona  amostras caso e controle
group<-c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
case_control30 <- pheno_data$Run[pheno_data$source_name == "Control_NPCs"&
                                     pheno_data$Day%in%group]
  
case_control30 <- append(case_control30, pheno_data$Run[pheno_data$source_name == paste0(name,"_NPCs")&
                                                            pheno_data$Day%in%group])
#cria matriz com logCPM e CPM
expLogCpm<- as.data.frame(logCPM[, colnames(logCPM)%in%case_control30])
expCpm<- as.data.frame(CPM[, colnames(CPM)%in%case_control30])
#remove os NA
expLogCpm<- na.omit(expLogCpm)
expCpm<- na.omit(expCpm)

# Estabelece os grupos do experimento de acordo com as amostras
targets <- as.factor(c(rep("control", 24),rep("lead30", 24)))

# renomeia as colunas da tabela de expressao de acordo com os grupos
colnames(expLogCpm) <- targets
colnames(expCpm) <- targets

# Obter a matriz design
lev <- levels(targets)
design <- model.matrix(~ 0 + targets)
colnames(design) <- c("control", "lead30")

# Modelo linear
fit <- lmFit(expLogCpm, design)
# Estabelecer os contrastes
# Neste momento, estabelecer as comparacoes entre os grupos 
# a fim de obter os valores de expressao diferencial
contrasts <- makeContrasts(lead30-control, levels=design)

# Estatistica bayesiana
ct.fit <- eBayes(contrasts.fit(fit, contrasts), trend = TRUE)
res.fit <- decideTests(ct.fit, method="global", adjust.method = "BH", p.value = 0.001)

# Combinar os resultados do limma em um unico dataset
SH.limma <- data.frame(logFC = ct.fit$coef, p.value = ct.fit$p.value, 
                       degenes = unclass(res.fit), stringsAsFactors = FALSE)


# Seleção dos genes diferencialmente expressos
features <- rowSums(res.fit!=0) > 0
features <- names(features)[features]

# Filtrar e obter a tabela dos genes diferencialmente expressos
DEexp <- expLogCpm[features, ]
DElimma <- SH.limma[rownames(DEexp), ]
DElimma <- DElimma[complete.cases(DElimma), ]

save(DEexp, DElimma, SH.limma, file = "DEexpTodos.RData")



# #seleciona lista de genes diferencialmente expressos em ambos os tempos
# DEGenes<-data.frame(ENSP=unique(c(transc[[1]]@DE$Protein,transc[[2]]@DE$Protein)))
# #transforma lista de proteinas em lista de genes
# DEGenes<-merge(DEGenes, dic_transcriptogramer, by = "ENSP" )
# DEGenes$ENSP<-NULL
# # exclui duplicidade
# DEGenes<-unique(DEGenes)

# #removendo genes não diferencialmente expressos da tabela CPM e logCPM
# expCpm<-expCpm[rownames(expCpm)%in%DEGenes$PROBE,]
# expLogCpm<-expLogCpm[rownames(expLogCpm)%in%DEGenes$PROBE,]

#removendo genes não diferencialmente expressos da tabela CPM e logCPM
expCpm<-expCpm[rownames(expCpm)%in%features,]
expLogCpm<-expLogCpm[rownames(expLogCpm)%in%features,]

save(expCpm, expLogCpm, file = "countsFiltered.RData")

ceaLogCpm <- cea(x= expLogCpm, sig=0.001, nper=10000, plot=T, p.adj.method = "BH")
ceaCpm <- cea(x= expCpm, sig=0.001, nper=10000, plot=T, p.adj.method = "BH")

save(ceaLogCpm,ceaCpm, file = "cea.RData")

l<-as.data.frame(ceaLogCpm)
c<-as.data.frame(ceaCpm)
for(i in 1:nrow(c)){
  for(j in 1:ncol(c)){
    if(abs(c[i,j] - l[i,j]) > 0.1){
      cat("Diferente em ",i," ",j,": ",
          c[i,j]," ",l[i,j],
          names(c)[i],names(l)[i],"\n")
    }
  }
}
