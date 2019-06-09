#Did you change it to your base location?
baseDir="~/LeadTest/"
setwd(baseDir)


GOdendogram <- function(object, cluster = 1, sizeIntervals = 10,
                        sizeMultiplier = 15, sizeBase = 100, cutoff = 0,
                        sizeLegend = TRUE, sizeLegendVertical = TRUE,
                        colorLegend = TRUE, colorLegendVertical = TRUE,
                        host = "127.0.0.1", port = 9091){
  terms <- object@Terms
  terms <- terms[terms$ClusterNumber %in% cluster,]
  if(nrow(terms)==0){
    stop("there are no enriched terms in this cluster!")
  }
  if(length(unique(terms$Annotated)) < sizeIntervals){
    message("** modifying sizeIntervals to avoid an error!")
    message("** the new value of sizeIntervals is ", length(unique(terms$Annotated)), "...")
    sizeIntervals <- length(unique(terms$Annotated))
  }
  GOs <- unique(terms$GO.ID)
  cutoff <- ceiling(nrow(terms)*cutoff)
  message("computing similarity matrix... step 1 of 4")
  mSim <- lapply(GOs, function(i, object, GOs){
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
  if(cutoff > 0){
    RS <- rowSums(mSim) - 1
    GOs <- names(sort(RS)[-1:-cutoff])
    mSim <- mSim[rownames(mSim) %in% GOs, colnames(mSim) %in% GOs]
    terms <- terms[terms$GO.ID %in% GOs,]
    GOs<-terms$GO.ID
  }
  message("generating graph from matrix... step 2 of 4")
  rtni <- RTN::tni.constructor(expData = mSim, regulatoryElements = GOs,
                               verbose = FALSE)
  rtni <- RTN::tni.permutation(rtni, verbose = FALSE)
  rtni <- RTN::tni.bootstrap(rtni, verbose = FALSE)
  rtni <- RTN::tni.dpi.filter(rtni, verbose = FALSE, eps = 0.05)
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
  pal <- RColorBrewer::brewer.pal(9, "YlOrRd")
  color_col <- grDevices::colorRampPalette(pal)(11)
  rate <- mapply(function(s, a){s/a}, terms$Significant, terms$Annotated)
#  rate <- rate/max(rate)
  rate <- rate/1
  igraph::V(graph)$rate <- NA
  igraph::V(graph)$rate[match(GOs,
                              igraph::V(graph)$name)] <- rate
  igraph::V(graph)$nodeAlias[match(GOs,
                                   igraph::V(graph)$name)] <- terms$Term
  graph <- RedeR::att.setv(g = graph, from = "rate", to = "nodeColor",
                           cols = color_col, na.col = "black", breaks = seq(0, 1, 0.1))
  rdp <- RedeR::RedPort(host = host, port = port)
  message("invoking RedeR... step 3 of 4")
  RedeR::calld(rdp)
  RedeR::resetd(rdp)
  myTheme <- list(zoom = 10)
  igraph::V(graph)$nodeFontSize[igraph::V(graph)$nodeFontSize == 20] <- 40
  suppressMessages(RedeR::addGraph(rdp, graph, theme = myTheme))
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
  message("relaxing nodes... step 4 of 4")
  
  teste<-data.frame(name=igraph::V(graph)$name,alias=igraph::V(graph)$nodeAlias,stringsAsFactors = F)
  #View(teste)
  
  RedeR::relax(rdp, p7 = 20, p8 = 20)
  message("done!")
  result<-list()
  result[[1]]<-rdp
  result[[2]]<-terms
  return(result)
  #return(rdp)
}

#rdp <- GOdendogram(t)
