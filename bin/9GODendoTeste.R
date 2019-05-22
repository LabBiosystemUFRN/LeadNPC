rm(list = ls())
setwd("/home/clovis/Dropbox/Chumbo/")
load("./transcriptograms/allTranscriptogramers80")
load("./Data/tfs.RData")
#load("./Data/colors.RData")
load("./Data/colorsSB.RData")

library(transcriptogramer)


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
# object <- t
 # clusters <- 7
 # nCores <- 1
 # sizeIntervals <- 10
 # sizeMultiplier <- 15
 # sizeBase <- 100

object<-transc[[1]]
clusters <- 11
nCores <- 4
sizeIntervals <- 10
sizeMultiplier <- 15
sizeBase <- 100
threshold = 0.5
GOdendogram <- function(object, clusters = 1, nCores = 1, sizeIntervals = 10,
                         sizeMultiplier = 15, sizeBase = 100, threshold = 0){
   terms <- object@Terms
   terms <- terms[terms$ClusterNumber %in% clusters,]
   GOs <- unique(terms$GO.ID)
   cl <- snow::makeSOCKcluster(nCores)
   on.exit(snow::stopCluster(cl))
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
   names(mSim) <- GOs
   mSim <- do.call("cbind", mSim)
   #mSim[mSim<=threshold]<-0
   message("generating graph from matrix...")
   rtni <- RTN::tni.constructor(expData = mSim, regulatoryElements = GOs,
                                verbose = T)
   rtni <- RTN::tni.permutation(rtni, verbose = T)
   rtni <- RTN::tni.bootstrap(rtni, verbose =T )
   rtni <- RTN::tni.dpi.filter(rtni, verbose = T)
   g <- RTN::tni.graph(rtni, tnet = "dpi", gtype = "amapDend", tfs = GOs)
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
   graph <- g$g
   igraph::V(graph)$nodeSize[match(GOs,
                                   igraph::V(graph)$name)] <- (sizeBase + (as.numeric(sizeRanks) * sizeMultiplier))
   #pal <- RColorBrewer::brewer.pal(9, "RdYlBu")
   pal <- RColorBrewer::brewer.pal(9, "YlOrRd")
   #color_col <- rev(grDevices::colorRampPalette(pal)(11))
   color_col <- grDevices::colorRampPalette(pal)(11)
   rate <- mapply(function(s, a){s/a}, terms$Significant, terms$Annotated)
   rate <- rate/max(rate)
   igraph::V(graph)$rate <- NA
   igraph::V(graph)$rate[match(GOs,
                               igraph::V(graph)$name)] <- rate
   igraph::V(graph)$nodeAlias[match(GOs,
                                    igraph::V(graph)$name)] <- terms$GO.ID
   # igraph::V(graph)$nodeAlias[match(GOs,
   #                                  igraph::V(graph)$name)] <- terms$Term
   graph <- RedeR::att.setv(g = graph, from = "rate", to = "nodeColor",
                            cols = color_col, na.col = "black", breaks = seq(0, 1, 0.1))
   rdp <- RedeR::RedPort()
   RedeR::calld(rdp)
   RedeR::resetd(rdp)
   myTheme <- list(zoom = 30)
   igraph::V(graph)$nodeFontSize[igraph::V(graph)$nodeFontSize == 20] <- 40
   RedeR::addGraph(rdp, graph, theme = myTheme)
   RedeR::addLegend.color(rdp, colvec = graph$legNodeColor$scale, size = 20,
                          labvec = graph$legNodeColor$legend, vertical = TRUE,
                          title = "Normalized term rate")
   RedeR::addLegend.size(rdp, sizevec = sort(unique(igraph::V(graph)$nodeSize)[-1]),
                         labvec = sizeDesc, vertical = TRUE,
                         title = paste("Terms size",nrow(terms)))
   #RedeR::relax(rdp, p1 = 50, p3 = 50, p8 = 20)
   RedeR::relax(rdp, p7 = 20, p8 = 20)
   return(rdp)
}

GOdendogram(object = object,clusters = 10,threshold = 0.9)

 object<-transc[[1]]
 clusters <- 10
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

 object@genesInTerm[[1]]
 