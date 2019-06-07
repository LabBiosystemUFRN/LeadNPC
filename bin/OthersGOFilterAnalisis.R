rm(list = ls())
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
setwd("/home/clovis/Dropbox/Chumbo/")

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
  #terms$X<-NULL
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
    #processa somente se nr de pontos for maior q 8
    # if(nrow(all)<=4){
    #   terms<-rbind(terms,allterms[(allterms$ClusterNumber==clNo),c(2,3,7,8)])
    # }else{
      #fit <- kmeans(all[,4:5], 4)
      #definindo o cluster + +
      #agr<-aggregate(all[,2:3],by=list(fit$cluster),FUN=mean)
      
      sdAnc<-sd(all$anc)
      mdAnc<-mean(all$anc)
      #GOs com mais filhos e menos pais
      all$score<-(1/(all$anc+1)*all$off)
      topGOs<-all[order(all$score,decreasing = T),][1:5,1]
      #qtScore<-quantile(all$score)
      #topGOs<-all$GO[all$score>qtScore[4]]
      #pega o cluster que possui em média mais filhos
      #clOk<-agr[agr[,3]==max(agr[,3]),1]
      #cluster mais equilbrado
      # agr<-apply(agr, MARGIN = 1,FUN = function(x){
      #     return(c(x[1],sd(c(x[2],x[3]))))
      #   })
      # agr<-t(agr)
      #clOk<-agr[agr[,2]==min(agr[,2]),1]
      #all <- data.frame(all,c= fit$cluster) 
      
      #cat("Cluster",clOk,"\n")
      
      # print(ggplot(data=all,aes(x=A,y=O))+
      #         labs(title=paste("Cluster", clNo,"- ok no ",clOk),
      #          x ="Pais", y = "Filhos")+
      #   geom_jitter(aes(color=factor(c))))
      print(ggplot()+
              labs(title=paste("Cluster", clNo),
               x ="Pais", y = "Filhos")+
        geom_jitter(data=all[all$GO%in%topGOs,],aes(x=anc,y=off),col="blue")+
        geom_point(data=all[!all$GO%in%topGOs,],aes(x=anc,y=off),col="red"))
        # terms<-rbind(terms,allterms[(allterms$ClusterNumber==clNo & 
      #                    allterms$GO.ID%in%all$GO[all$c == clOk]),c(2,3,7,8)])
      terms<-rbind(terms,allterms[(allterms$ClusterNumber==clNo & 
                                     allterms$GO.ID%in%topGOs),c(1,2,6,7)])
  #  }
    
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





# #*************** código Diego
# GOTreeRedeR <- function(object, clusters = 1, ontology = "BP"){
#   data <- transcriptogramer::Terms(object)
#   data <- data[data$ClusterNumber %in% clusters, ]
#   rootName <- hashTable <- rootGO <- NULL
#   if(ontology == "BP"){
#     rootName <- "biological process"
#     rootGO <- "GO:0008150"
#     hashTable <- as.list(GO.db::GOBPPARENTS)
#   }else if(ontology == "CC"){
#     rootName <- "cellular component"
#     rootGO <- "GO:0005575"
#     hashTable <- as.list(GO.db::GOCCPARENTS)
#   }else if(ontology == "MF"){
#     rootName <- "molecular function"
#     rootGO <- "GO:0003674"
#     hashTable <- as.list(GO.db::GOMFPARENTS)
#   }
#   df <- data.frame()
#   GOs <- unique(data$GO.ID)
#   invisible(sapply(GOs, function(i) {
#     if(i==rootGO){
#       next
#     }
#     nParents <- 1
#     parents <- hashTable[[which(names(hashTable)==i)]]
#     if(any(parents %in% GOs)){
#       idx <- which(parents %in% GOs)
#       if(length(idx)>1){
#         # More than one parent
#         parents <- parents[idx]
#         nParents <- length(parents)
#         parents <- sapply(parents, function(term){GO.db::GOTERM[[term]]@Term})
#       }else{
#         parents <- parents[idx]
#         parents <- GO.db::GOTERM[[parents]]@Term
#       }
#     }else{
#       parents <- rootName
#     }
#     aux <- data.frame(parent = parents, child = rep(GO.db::GOTERM[[i]]@Term, nParents),
#                       ID = rep(i, nParents),
#                       stringsAsFactors = F)
#     df <<- rbind(df, aux)
#     return(NULL)
#   }))
#   return(df)
#   rdp <- RedeR::RedPort()
#   RedeR::calld(rdp)
#   possibleError<-tryCatch(g <- igraph::graph.data.frame(df),
#            error=function(e) e)
#   if(inherits(possibleError, "error")) {
#     cat("Erro!")
#     return(rdp)
#   }
#   #tamanho do nó
#   # igraph::V(g)$nodeSize <- vapply(igraph::V(g)$name, function(i){
#   #   return(as.integer(igraph::neighborhood.size(g, igraph::vcount(g), i, "out")+10))
#   # }, integer(1))
#   igraph::E(g)$edgeColor <- "grey80"
#   igraph::V(g)$nodeLineColor <- "grey80"
#   g <- RedeR::att.mapv(g, df, 2)
# #  igraph::V(g)$nodeAlias <- igraph::V(g)$name
#   igraph::V(g)$nodeAlias <- igraph::V(g)$ID
#   igraph::V(g)$pValue <- data[match(igraph::V(g)$ID, data$GO.ID), 6]
#   aux <- na.omit(igraph::V(g)$pValue)
#   igraph::V(g)$pValue <- as.factor(igraph::V(g)$pValue)
#   pal <- RColorBrewer::brewer.pal(9, "RdYlBu")
#   color_col <- grDevices::colorRampPalette(pal)(length(unique(data$pValue)))
#   g <- RedeR::att.setv(g = g, from = "pValue", to = "nodeColor", cols = color_col,
#                        na.col = "grey80", breaks = seq(1, length(color_col)))
#   #RedeR::addGraph(rdp, g)
#   #RedeR::addGraph(rdp, g, igraph::layout_as_tree(g,flip.y = F))
#   #RedeR::addLegend.color(rdp, colvec = g$legNodeColor$scale, size = 15,
#   #                       labvec = sort(unique(aux)),
#   #                       title = "Adjusted p-values")
#   #RedeR::relax(rdp)
#   #return(rdp)
#   return(g)
# }
# library(igraph)
# load("./transcriptograms/allTranscriptogramers80")
# 
# i =1
# #for(i in 1:2){
#   # #qtd de clusteres
#   # clusters<-unique(allterms$ClusterNumber)
#   
#   #object <- transc[[i]]
#   #t=c(2,1)
# clusteres<-data.frame(T1=c(1,	2,	3,	4,	5,	6,	7,	8,	9,	10,	
#                            11,	12,	13,	14,	15,	16,	17,NA,NA),
#                       T2=c(2,	3,	3,	4,	5,	6,	9,	10,	11,	11,	11,	11,	11,	11,	12,	13,	13,7,8))
# save(clusteres,file = "./Data/clusteres.RData")
#   t2 = 13
#   rdp <- RedeR::RedPort()
#   RedeR::calld(rdp)
#   for(t2 in 7:13){
#     if(t2%in%c(7,8)){
#       df2 <- GOTreeRedeR(transc[[2]],t2)
#       df2$itv<-2
#       df3<-df2
#     }else  if(t2==13){
#       df1 <- GOTreeRedeR(transc[[1]],t1)
#       df1$itv<-1
#       df3<-df1
#     }else{
#       t1<-clusteres$T1[clusteres$T2==t2]
#       #RedeR::resetd(rdp)
#       #rdp <- GOTreeRedeR(object,t)
#       df1 <- GOTreeRedeR(transc[[1]],t1)
#       df1$interval<-1
#       df2 <- GOTreeRedeR(transc[[2]],t2)
#       df2$interval<-2
#       df3<-merge(df1,df2,by=c("parent","child","ID"),all=T)
#       df3$interval.x[is.na(df3$interval.x)]<-0
#       df3$interval.y[is.na(df3$interval.y)]<-0
#       df3$itv<-df3$interval.x+df3$interval.y
#       df3$interval.x<-NULL
#       df3$interval.y<-NULL
#     }
#     g3 <- igraph::graph.data.frame(df3)
#     g3 <- RedeR::att.mapv(g3, df3, 2)
#     V(g3)$itv
#     # ng1<-V(g1)$ID
#     # ng2<-V(g2)$ID
#     V(g3)$nodeColor[V(g3)$itv == 1]<-"#00FF00"
#     V(g3)$nodeColor[V(g3)$itv == 2]<-"#0000FF"
#     V(g3)$nodeColor[V(g3)$itv == 3]<-"#FF00FF"
#     igraph::V(g3)$nodeAlias <- igraph::V(g3)$ID
#     V(g3)$nodeSize<-10
#     RedeR::resetd(rdp)
#     RedeR::addGraph(rdp, g3, igraph::layout_as_tree(g3,flip.y = F))
#     readline(prompt=paste0("Cl",t2, ". Enter continua... "))
#   }
# #}
# 
# 

