rm(list = ls())
setwd("/home/clovis/Dropbox/Chumbo/")
load("./transcriptograms/allTranscriptogramers80")
load("./Data/tfs.RData")
#load("./Data/colors.RData")
load("./Data/colorsSB.RData")

library(transcriptogramer)
source(file = "./bin/GOdendogram.R")

# t <- transcriptogramPreprocess(association = association, ordering = Hs900)
# t <- transcriptogramStep1(object = t, expression = GSE9988,
#                           dictionary = GPL570, nCores = T)
# radius(object = t) <- 80
# t <- transcriptogramStep2(object = t, nCores = T)
# levels <- c(rep(FALSE, 3), rep(TRUE, 3))
# t <- differentiallyExpressed(object = t, levels = levels, pValue = 0.01)
# t <- clusterEnrichment(object = t, species = HsBPTerms,
#                        pValue = 0.005, nCores = T)
# #save(t, file = "transcriptogram.RData", compress = "xz")
# genesInTerm <- t@genesInTerm
# plot <- enrichmentPlot(t, T)
# terms <- Terms(t)
# length(unique(terms$GO.ID))
# 
# load("transcriptogram.RData")
# 
tempo=1
object<-transc[[tempo]]
clusters <- 6
nCores <- 10
sizeIntervals <- 10
sizeMultiplier <- 15
sizeBase <- 100

 # GOdendogram <- function(object, clusters = 1, nCores = 1, sizeIntervals = 10,
 #                         sizeMultiplier = 15, sizeBase = 100){
 #   #corte em 10% da janela
 #   terms <- object@Terms
 #   terms <- terms[terms$ClusterNumber %in% clusters & terms$Significant>=16,]
 #   if(nrow(terms)==0){
 #     cat("Nada a fazer...")
 #     return(0)
 #   }
 #   GOs <- unique(terms$GO.ID)
 #   cl <- snow::makeSOCKcluster(nCores)
 #   on.exit(snow::stopCluster(cl))
 #   message("computing similarity matrix...")
 #   mSim <- snow::parLapply(cl, GOs, function(i, object, GOs){
 #     v <- sapply(GOs, function(j){
 #       union <- length(unique(c(object@genesInTerm[[i]], object@genesInTerm[[j]])))
 #       if(union==0){
 #         return(0)
 #       }
 #       inter <- sum(object@genesInTerm[[i]] %in% object@genesInTerm[[j]])
 #       jaccard <- inter/union
 #       return(jaccard)
 #     })
 #     return(v)
 #   }, object, GOs)
 #   names(mSim) <- GOs
 #   mSim <- do.call("cbind", mSim)
 #   message("generating graph from matrix...")
 #   rtni <- RTN::tni.constructor(expData = mSim, regulatoryElements = GOs,
 #                                verbose = FALSE)
 #   rtni <- RTN::tni.permutation(rtni, verbose = FALSE)
 #   rtni <- RTN::tni.bootstrap(rtni, verbose = FALSE)
 #   rtni <- RTN::tni.dpi.filter(rtni, verbose = FALSE)
 #   g <- RTN::tni.graph(rtni, tnet = "dpi", gtype = "amapDend", tfs = GOs)
 #   size <- terms$Annotated
 #   size <- classInt::classIntervals(size, sizeIntervals)
 #   aux <- as.character(size$brks)
 #   sizeDesc <- vapply(1:(length(size$brks) - 1), function(i){
 #     if(i == (length(size$brks) - 1)){
 #       return(paste0("[", aux[i], ", ", aux[i+1], "]"))
 #     }
 #     return(paste0("[", aux[i], ", ", aux[i+1], ")"))
 #   }, character(1))
 #   rm(aux)
 #   sizeRanks <- cut(size$var, breaks = size$brks, labels = as.character(1:sizeIntervals),
 #                    include.lowest = TRUE)
 #   size <- size$var
 #   graph <- g$g
 #   igraph::V(graph)$nodeSize[match(GOs,
 #                                   igraph::V(graph)$name)] <- (sizeBase + (as.numeric(sizeRanks) * sizeMultiplier))
 #   #pal <- RColorBrewer::brewer.pal(9, "RdYlBu")
 #   pal <- RColorBrewer::brewer.pal(9, "YlOrRd")
 #   #color_col <- rev(grDevices::colorRampPalette(pal)(11))
 #   color_col <- grDevices::colorRampPalette(pal)(11)
 #   rate <- mapply(function(s, a){s/a}, terms$Significant, terms$Annotated)
 #   rate <- rate/max(rate)
 #   igraph::V(graph)$rate <- NA
 #   igraph::V(graph)$rate[match(GOs,
 #                               igraph::V(graph)$name)] <- rate
 #   igraph::V(graph)$nodeAlias[match(GOs,
 #                                    igraph::V(graph)$name)] <- terms$Term
 #   graph <- RedeR::att.setv(g = graph, from = "rate", to = "nodeColor",
 #                            cols = color_col, na.col = "black", breaks = seq(0, 1, 0.1))
 #   rdp <- RedeR::RedPort()
 #   RedeR::calld(rdp)
 #   RedeR::resetd(rdp)
 #   myTheme <- list(zoom = 3)
 #   igraph::V(graph)$nodeFontSize[igraph::V(graph)$nodeFontSize == 20] <- 40
 #   RedeR::addGraph(rdp, graph, theme = myTheme)
 #   RedeR::addLegend.color(rdp, colvec = graph$legNodeColor$scale, size = 20,
 #                          labvec = graph$legNodeColor$legend, vertical = TRUE,
 #                          title = "Normalized term rate")
 #   RedeR::addLegend.size(rdp, sizevec = sort(unique(igraph::V(graph)$nodeSize)[-1]),
 #                         labvec = sizeDesc, vertical = TRUE,
 #                         title = paste("Terms size",nrow(terms)))
 #   #RedeR::relax(rdp, p1 = 50, p3 = 50, p8 = 20)
 #   #RedeR::relax(rdp, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   RedeR::relax(rdp, p1=220,p3=2000, p7 = 20, p8 = 5000)
 #   Sys.sleep(5)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   RedeR::relax(rdp, p1=20, p7 = 20, p8 = 20)
 #   
 #   return(rdp)
 # }

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
 source(file = "./bin/GOdendogramErro.R.R")
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
#View(terms)
 #cl7
 #descTerms<-c("negative regulation of nucleic acid-temp...","negative regulation of RNA biosynthetic ...","positive regulation of cellular process","negative regulation of metabolic process","negative regulation of cell cycle","negative regulation of gene expression, ...","gene silencing","regulation of cyclin-dependent protein s...","positive regulation of RNA metabolic pro...","positive regulation of cellular biosynth...","endoplasmic reticulum to cytosol transpo...","regulation of gene silencing","chromatin organization","protein modification by small protein co...","cell cycle process","cellular response to stress","cellular macromolecule metabolic process","protein catabolic process","chromosome segregation","proteolysis involved in cellular protein...")
 #cl8
 #descTerms<-c("ncRNA transcription","snRNA transcription by RNA polymerase II","nucleobase-containing small molecule met...","nucleoside phosphate biosynthetic proces...","ribose phosphate metabolic process","cellular aromatic compound metabolic pro...","snRNA metabolic process","transcription by RNA polymerase III","ncRNA metabolic process")
 #cl9
 #descTerms<-c("nuclear transport","mRNA metabolic process","mRNA processing","RNA processing","protein folding","protein import into nucleus","DNA-templated transcription, termination","RNA localization","RNA splicing","mRNA export from nucleus")
 #cl10
 #descTerms<-c("mRNA metabolic process","peptide metabolic process","ribonucleoprotein complex biogenesis","translation","translational initiation","ncRNA metabolic process","RNA catabolic process","mRNA catabolic process","posttranscriptional regulation of gene e...")
 
 # duplicated(descTerms)
 # length(descTerms)
 # nodeList<-terms$GO.ID[terms$Term%in%descTerms]
 # #nodeList<-terms$Term[terms$Term%in%descTerms]
 # #descTerms[!descTerms%in%nodeList]
 # View(terms[terms$Term%in%descTerms,])
 # 
 # descTerms[order(descTerms)]
 

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

