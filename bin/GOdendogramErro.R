# # load()
# tempo=1
# cluster =7
# 
# object<-transc[[tempo]]
# #object <- t
# #cluster <- 1
# nCores <- 10
# sizeIntervals <- 10
# sizeMultiplier <- 15
# sizeBase <- 100
# sizeLegend <- TRUE
# colorLegend <- TRUE
# sizeLegendVertical <- TRUE
# colorLegendVertical <- TRUE
# cutoff <- 0.25

GOdendogram <- function(object, cluster = 1, nCores = 1, sizeIntervals = 10,
                        sizeMultiplier = 15, sizeBase = 100, cutoff = 0,
                        sizeLegend = TRUE, sizeLegendVertical = TRUE,
                        colorLegend = TRUE, colorLegendVertical = TRUE, maxRate=NA){
  terms <- object@Terms
  terms <- terms[terms$ClusterNumber %in% cluster,]
  GOs <- unique(terms$GO.ID)
  cutoff <- ceiling(nrow(terms)*cutoff)
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
  if(cutoff > 0){
    RS <- rowSums(mSim) - 1
    GOs <- names(sort(RS)[-1:-cutoff])
    mSim <- mSim[rownames(mSim) %in% GOs, colnames(mSim) %in% GOs]
    terms <- terms[terms$GO.ID %in% GOs,]
  }
  message("generating graph from matrix...")
  rtni <- RTN::tni.constructor(expData = mSim, regulatoryElements = GOs,
                               verbose = FALSE)
  rtni <- RTN::tni.permutation(rtni, verbose = FALSE)
  rtni <- RTN::tni.bootstrap(rtni, verbose = FALSE)
  rtni <- RTN::tni.dpi.filter(rtni, verbose = FALSE, eps = 0.05)
  #g <- RTN::tni.graph(rtni, tnet = "dpi", gtype = "amapDend", tfs = GOs)
  grafo = function(rtni){
    g<-RTN::tni.graph(rtni, tnet = "dpi", gtype = "amapDend", tfs = GOs)
    graph <- g$g
    return(graph)
  }
  graph<-tryCatch( grafo(rtni) ,
               error = function(e){
                 graph<-graph_from_adjacency_matrix(rtni@results$tn.dpi, mode = c("undirected"))
                 cat("Warning: no conections detectec\n")
                 return(graph)
               })
  
  
  size <- terms$Annotated
  if(length(size) < sizeIntervals){
    sizeIntervals = length(size)
  }
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
  igraph::V(graph)$nodeSize[match(GOs,
                                  igraph::V(graph)$name)] <- (sizeBase + (as.numeric(sizeRanks) * sizeMultiplier))
  #pal <- RColorBrewer::brewer.pal(9, "RdYlBu")
  pal <- RColorBrewer::brewer.pal(9, "YlOrRd")
  #color_col <- rev(grDevices::colorRampPalette(pal)(11))
  color_col <- grDevices::colorRampPalette(pal)(11)
  rate <- mapply(function(s, a){s/a}, terms$Significant, terms$Annotated)
  if(is.na(maxRate)){
    rate <- rate/max(rate)
  }else{
    #cria mesma escala de cores para a figura 2
    rate <- rate/maxRate
  }
  igraph::V(graph)$rate <- NA
  igraph::V(graph)$rate[match(GOs,
                              igraph::V(graph)$name)] <- rate
  igraph::V(graph)$nodeAlias[match(GOs,
                                   igraph::V(graph)$name)] <- terms$Term
  graph <- RedeR::att.setv(g = graph, from = "rate", to = "nodeColor",
                           cols = color_col, na.col = "black", breaks = seq(0, 1, 0.1))
  rdp <- RedeR::RedPort()
  RedeR::calld(rdp)
  RedeR::resetd(rdp)
  myTheme <- list(zoom = 10)
  if(is.null(igraph::V(graph)$nodeFontSize)){
    igraph::V(graph)$nodeFontSize <- 40
  }
  igraph::V(graph)$nodeFontSize[igraph::V(graph)$nodeFontSize == 20] <- 40
  RedeR::addGraph(rdp, graph, theme = myTheme)
  if(colorLegend){
    RedeR::addLegend.color(rdp, colvec = graph$legNodeColor$scale, size = 20,
                           labvec = graph$legNodeColor$legend, vertical = colorLegendVertical,
                           title = "Normalized term rate")
  }
  if(sizeLegend){
    RedeR::addLegend.size(rdp, sizevec = sort(unique(igraph::V(graph)$nodeSize)[-1]),
                          labvec = sizeDesc, vertical = sizeLegendVertical,
                          title = paste0("Term size (", nrow(terms), " terms)"))
  }
  #RedeR::relax(rdp, p1 = 50, p3 = 50, p8 = 20)
  
  RedeR::relax(rdp, p7 = 20, p8 = 20)
  
  return(rdp)
}

#rdp <- GOdendogram(t, cutoff = 0.25, colorLegendVertical = FALSE, nCores = 4)
