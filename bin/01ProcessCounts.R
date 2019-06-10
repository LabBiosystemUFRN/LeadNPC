rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadNPC/"
setwd(baseDir)
source(file = "bin/00base.R")

library("edgeR")
library(biomaRt)

setwd(paste0(baseDir,"Data/"))

samples<-read.table("allCountsHisat.txt",stringsAsFactors = F,header = T, sep = "\t")
x<-samples[,7:59]
#Cria dicionario para transcriptogramer ----
#acerta o nome das coluna de amostra
dic_transcriptogramer<-data.frame(PROBE=samples$Geneid)

#Ensemble gene ID passa a ser o probe associado ao peptide id
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

#lista<-listAttributes(mart)
lista <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), mart = mart)

dic_transcriptogramer<-merge(dic_transcriptogramer,lista, 
                             by.x="PROBE",by.y = "ensembl_gene_id")

dic_transcriptogramer<-subset(dic_transcriptogramer, ensembl_peptide_id != "")
dic_transcriptogramer<-subset(dic_transcriptogramer, ensembl_peptide_id != " ")
dic_transcriptogramer<-data.frame(ENSP=as.character(dic_transcriptogramer$ensembl_peptide_id),
                                  PROBE = as.character(dic_transcriptogramer$PROBE),stringsAsFactors = F)
colnames(dic_transcriptogramer)<-c("ENSP","PROBE")

#acrescenta o prefixo 9606.
dic_transcriptogramer$ENSP<-gsub("ENSP","9606.ENSP",
                                 dic_transcriptogramer$ENSP,fixed = T)

load(file =  "phenodata_Orig.RData")

# processamento das amostras ----
#Corte de amostras com menos de 1 cpm
#ref: RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR
cpm <- cpm(samples[,7:59])
keep.exprs <-rowSums(cpm>1)>=3
sum(keep.exprs)
samples <- samples[keep.exprs,]

#mantém no nome das colunas de amostras somente o ensbl gene ID
colnames(samples)<-c(colnames(samples[,1:6]),substr(colnames(samples[,7:59]),1,10))

#Converte para logCPM
niente<-as.data.frame(cpm(samples[,c(7:59)],log = T))
rownames(niente)<-samples$Geneid

logCPM<-as.matrix(niente)

niente<-as.data.frame(cpm(samples[,c(7:59)],log = F))
rownames(niente)<-samples$Geneid
CPM<-as.matrix(niente)

save(dic_transcriptogramer,logCPM,pheno_data, file = "counts.RData")
save(CPM,file = "CPM.RData")
# suplemento ----
#plota diferença entre dataset filtrado e o sem filtro
library(RColorBrewer)
lcpm <- cpm(x, log=T)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.6),  las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

cpm <- cpm(x, log=TRUE)
keep.exprs <-rowSums(cpm>1)>=3
x1 <- x[keep.exprs,]
lcpm <- cpm(x1, log=TRUE)

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.6), las=2,
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2 : nsamples){
  den <- density(lcpm[,i])
  lines(den$x , den$y, col=col[i], lwd=2)
}

