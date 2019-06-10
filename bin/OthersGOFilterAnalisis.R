rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadNPC/"
setwd(baseDir)
source(file = "bin/00base.R")

library(topGO)
library(tibble)
library(stats)
library(ggplot2)


# consulta<-c("GO:0002228","GO:0002715","GO:0006497","GO:0006505","GO:0006506",
#             "GO:0006643","GO:0008544","GO:0008654","GO:0009247","GO:0009913",
#             "GO:0015695","GO:0016254","GO:0016255","GO:0019538","GO:0031424",
#             "GO:0042157","GO:0042158","GO:0045017","GO:0046467","GO:0046486",
#             "GO:0050954","GO:0050957","GO:0055085","GO:0070268","GO:0090407",
#             "GO:1901137","GO:1903509")
 consulta<-c("GO:0051606","GO:0009593","GO:0006629","GO:0007186","GO:0008610","GO:0050906","GO:0044255","GO:0044281","GO:0003008","GO:0006082","GO:0007606","GO:0055114","GO:0042221","GO:1901615","GO:0044283","GO:0008202","GO:0006694","GO:1901617","GO:0050877","GO:0016054","GO:0050907","GO:0044242","GO:0007600","GO:0006637","GO:0097164","GO:0016053","GO:0035383","GO:0032501","GO:0050911","GO:0007608","GO:0006672","GO:0016042","GO:0006790","GO:0044282","GO:0006665","GO:0042430","GO:0046513","GO:0009072","GO:0042537","GO:0044106","GO:0030148","GO:0034440","GO:0046890","GO:0042436","GO:0009063","GO:0051186","GO:0009074","GO:0017144")

#ancestrais contidos em GO.db
ancT <- as.list(GOBPANCESTOR)
#remove NA
ancT[is.na(ancT)]<-"Nada"

#descendentes contidos em GO.db
offT <- as.list(GOBPOFFSPRING)
#remove NA, mas mantém item na lista
offT[is.na(offT)]<-"Nada"
intervalo=1

for(intervalo in 1:2){
  #carrega arquivo de enrich GO
  allterms<-read.table(paste0("/home/clovis/Dropbox/Chumbo/terms/AllW80PCAGr",intervalo,"Lead30.csv"),sep = "\t",stringsAsFactors = F)
  #remove termos duplicados dentro dos clusteres
  goDup<-allterms$GO.ID[duplicated(allterms$GO.ID)]
  allterms<-allterms[!(allterms$GO.ID%in%goDup),]
  #qtd de clusteres
  clusters<-unique(allterms$ClusterNumber)
  
  #prepara dataframe de resultados
  terms<- allterms[1,]
  terms<-terms[-1,]
  terms$Annotated<-NULL
  terms$Significant<-NULL
  terms$Expected<-NULL
  terCol<-colnames(terms)
  
  clNo=6
  for(clNo in clusters){
    #extrai GOs contidas nos clusteres
    consulta<-allterms$GO.ID[allterms$ClusterNumber == clNo]
    
    #filtra por consulta
    anc<-ancT[names(ancT)%in%consulta]
    off<-offT[names(offT)%in%consulta]
    
    #filtra lista de ancestrais
    anc<-lapply(anc, function(x){
      x[x%in%consulta]
      })
    #filtra lista de descendentes
    off<-lapply(off, function(x){
      x[x%in%consulta]
    })
    
    #conta itens
    ancestor<-as.data.frame(sapply(X = anc, FUN = function(x){
      return(length(x))
    }))
    
    offSpring<-as.data.frame(sapply(X = off, FUN = function(x){
      return(length(x))
    }))
    ancestor<-tibble::rownames_to_column(ancestor)
    offSpring<-tibble::rownames_to_column(offSpring)
    
    colnames(ancestor)<-c("GO","Anc")
    colnames(offSpring)<-c("GO","Off")
    
    all<-merge(ancestor,offSpring,all=T, by="GO")
    all$Off[is.na(all$Off)]<-0
    all$Anc[is.na(all$Anc)]<-0
    
    all<-cbind(all,scale(all[,2:3]))
    colnames(all)<-c("GO","anc","off","Anc","Off")
    all$A<-all$anc-median(all$anc)
    all$O<-all$off-median(all$off)
      
      sdAnc<-sd(all$anc)
      mdAnc<-mean(all$anc)
      #GOs com mais filhos e menos pais
      all$score<-(1/(all$anc+1)*all$off)
      topGOs<-all[order(all$score,decreasing = T),][1:5,1]
      print(ggplot()+
              labs(title=paste("Cluster", clNo),
               x ="Pais", y = "Filhos")+
        geom_jitter(data=all[all$GO%in%topGOs,],aes(x=anc,y=off),col="blue")+
        geom_point(data=all[!all$GO%in%topGOs,],aes(x=anc,y=off),col="red"))
      terms<-rbind(terms,allterms[(allterms$ClusterNumber==clNo & 
                                     allterms$GO.ID%in%topGOs),c(1,2,6,7)])
    
  }
  
  intTerms<-c("calcium","zinc")
  for(cont in intTerms){
    terms<-rbind(terms,allterms[grep(cont,allterms$Term),c(1,2,6,7)])
  }
  terms<-terms[order(terms$ClusterNumber),]
  write.table(terms, file = paste0("./terms/topTermsInterval",intervalo,".csv"), 
              sep="\t",
              row.names = F)
  
}





