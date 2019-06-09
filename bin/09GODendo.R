rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadTest/"
setwd(baseDir)
source(file = "bin/00base.R")

load("./Data/allTranscriptogramers80")
load("./Data/tfs.RData")
load("./Data/colors.RData")

library(transcriptogramer)
source(file = "./bin/GOdendogram.R")

 
tempo=1
object<-transc[[tempo]]
clusters <- 6
nCores <- 10
sizeIntervals <- 10
sizeMultiplier <- 15
sizeBase <- 100


 clusters <- 6
 nCores <- 1
 sizeIntervals <- 10
 sizeMultiplier <- 15
 sizeBase <- 100
 
 GOfilter <- function(object, clusters = 1, nCores = 1, sizeIntervals = 10,
                         sizeMultiplier = 15, sizeBase = 100){
   terms <- object@Terms
   terms <- terms[terms$ClusterNumber %in% clusters,]
   GOs <- unique(terms$GO.ID)
   cl <- snow::makeSOCKcluster(nCores)
   on.exit(snow::stopCluster(cl))
   
   message(paste("computing similarity matrix cluster ",clusters," ..."))
   mSim <- snow::parLapply(cl, GOs, function(i, object, GOs){
     v <- sapply(GOs, function(j){
       union <- length(unique(c(object@genesInTerm[[i]], object@genesInTerm[[j]])))
       if(union==0){
         return(0)
       }
       inter <- sum(object@genesInTerm[[i]] %in% object@genesInTerm[[j]])
       jaccard <- inter/union
       return(jaccard)
     })
     return(v)
   }, object, GOs)
   names(mSim) <- GOs
   mSim <- do.call("cbind", mSim)
   message(paste("generating graph from matrix cluster ",clusters," ..."))
   if(nrow(terms)>=20){
     rtni <- RTN::tni.constructor(expData = mSim, regulatoryElements = GOs,
                                  verbose = FALSE)
     rtni <- RTN::tni.permutation(rtni, verbose = FALSE)
     rtni <- RTN::tni.bootstrap(rtni, verbose = FALSE)
     rtni <- RTN::tni.dpi.filter(rtni, verbose = FALSE)
     #g <- RTN::tni.graph(rtni, tnet = "dpi", gtype = "amapDend", tfs = GOs)
     size <- terms$Annotated
     size <- classInt::classIntervals(size, sizeIntervals)
     aux <- as.character(size$brks)
     sizeDesc <- vapply(1:(length(size$brks) - 1), function(i){
       if(i == (length(size$brks) - 1)){
         return(paste0("[", aux[i], ", ", aux[i+1], "]"))
       }
       return(paste0("[", aux[i], ", ", aux[i+1], ")"))
     }, character(1))
     rm(aux)
     sizeRanks <- cut(size$var, breaks = size$brks, labels = as.character(1:sizeIntervals),
                      include.lowest = TRUE)
     size <- size$var
     rate <- mapply(function(s, a){s/a}, terms$Significant, terms$Annotated)
     rate <- rate/max(rate)
     terms$rate<-rate
   }else{
     terms$rate<-rep(NA,nrow(terms))
   }
   return(terms)
 }
 
 tempo=1
 object<-transc[[tempo]]
 
 #determina maxRate
 clusteres<-c(5,6,11)
 mRate<-c()
 cutoff=0.25
 i=2
 cluster=7
 nCores = 10
 for(i in clusteres){
   terms <- object@Terms
   terms <- terms[terms$ClusterNumber %in% cluster,]
   GOs <- unique(terms$GO.ID)
   cutoff2 <- ceiling(nrow(terms)*cutoff)
   cl <- snow::makeSOCKcluster(nCores)
   #on.exit(snow::stopCluster(cl))
   message("computing similarity matrix...")
   mSim <- snow::parLapply(cl, GOs, function(i, object, GOs){
     v <- sapply(GOs, function(j){
       union <- length(unique(c(object@genesInTerm[[i]], object@genesInTerm[[j]])))
       if(union==0){
         return(0)
       }
       inter <- sum(object@genesInTerm[[i]] %in% object@genesInTerm[[j]])
       jaccard <- inter/union
       return(jaccard)
     })
     return(v)
   }, object, GOs)
   snow::stopCluster(cl)
   names(mSim) <- GOs
   mSim <- do.call("cbind", mSim)
   if(cutoff2 > 0){
     RS <- rowSums(mSim) - 1
     GOs <- names(sort(RS)[-1:-cutoff2])
     mSim <- mSim[rownames(mSim) %in% GOs, colnames(mSim) %in% GOs]
     terms <- terms[terms$GO.ID %in% GOs,]
   }
   
   rate <- mapply(function(s, a){s/a}, terms$Significant, terms$Annotated)
   mRate <- c(mRate,max(rate))
 }
 
 maxRate<-max(mRate)
 
 tempo=1
 object<-transc[[tempo]]
 
 cluster = 5
 #maxRate=1
 #source(file = "./bin/GOdendogramErro.R.R")
 result <- tryCatch(GOdendogram(object, 
                               cluster = cluster, 
                               #nCores = 10,
                               cutoff = 0, 
                               colorLegendVertical = FALSE),
                             #maxRate = maxRate),
                   error = function(e) e)
 terms<-result[[2]]
 rdp<-result[[1]]
 
 RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000,ps=T)
 #readline(prompt=paste0("t",tempo, "c",clust," Acerte a topologia. Enter continua... "))
 RedeR::relax(rdp, p1=100, p7 = 20, p8 = 20,ps=T)
 RedeR::relax(rdp, p1=30, p7 = 20, p8 = 20,ps=T)

 RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000,ps=T)
 #readline(prompt=paste0("t",tempo, "c",clust," Acerte a topologia. Enter continua... "))
 RedeR::relax(rdp, p1=100, p7 = 20, p8 = 20,ps=T)
 RedeR::relax(rdp, p1=30, p7 = 20, p8 = 20,ps=T)
 
 
 RedeR::selectNodes(rdp,nodeList)
 
library(RedeR)
 rdp <- RedPort('MyPort') 
 calld(rdp)

 #loop dendo ----
 for(tempo in 1:2){
   object<-transc[[tempo]]
   maxCl<-length(object@clusters)
   for (clust in 1:maxCl) {
     rdp <- tryCatch(GOdendogram(object, 
                                 cluster = clust, 
                                 nCores = 10,
                                 cutoff = 0.25, 
                                 colorLegendVertical = FALSE),
                     error = function(e) e)
     readline(prompt=paste0("t",tempo, "c",clust,". Enter continua... "))
   }
 }
 
 
  
 
 # "GO:0010257" "GO:0032981"
 
 
for(tempo in 1:2){
  for (clust in 1:11) {
    object<-transc[[tempo]]
    if(clust ==1){
      terms <- GOfilter(object, clusters = clust)
    }else{
      terms <- rbind(terms,GOfilter(object, clusters = clust))
    }
  }
  write.csv(terms, 
            file = paste0("/home/clovis/Dropbox/Chumbo/terms/termsInterval",
                          tempo,".csv"))
}

 
 
 library(transcriptogramer)
 object<-transc[[1]]
GOdendogram(object,cutoff=0.25,cluster = 11) 
RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000,ps=T)
RedeR::relax(rdp, p1=100, p7 = 20, p8 = 20,ps=T)
RedeR::relax(rdp, p1=30, p7 = 20, p8 = 20,ps=T)

