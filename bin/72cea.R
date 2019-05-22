rm(list = ls())

# COEXPRESSION ANALYSIS ---------------------------------------------------
library(RedeR)
library(igraph)
setwd("/home/clovis/Dropbox/Chumbo/Data/")


# Load arquivos
load("clusters12.RData")
#load("exp_diff_3GRUPOS.RData")
load("./allTranscriptogramers80")
load(file = "./counts.RData")
load("tfs.RData")
load("cea.RData")

DElimma2<-rbind(transc[[1]]@DE,transc[[2]]@DE)
DElimma2<-transc[[2]]@DE
DElimma2<-merge(DElimma2,dic_transcriptogramer,by.x = "Protein", by.y="ENSP")

clusters <- clusters12

# Abrir porta do RedeR
rdp <- RedPort()
calld(rdp)

# ESCOLHER OS CLUSTERS FORMADOS A PARTIR DO GRUPO 3 -----------------------

## PLOTAR LOGFC ------------------------------------------------------------

## Plotar as redes apenas com os valores de expressao que relativos a cada um
## dos grupos

cluster=3
logfc_grupo="logFC.Lead30_5...Control_5"
cut_level = 2
breaks = seq(-5, 5, 0.5)

plot_nests <- function(cluster, logfc_grupo, cut_level = 2, breaks = seq(-5, 5, 0.5)){
  print(cluster)
  
  print(paste0("Rede - cluster ", cluster))
  
  #Filtra porteinas e probes por clusters
  proteins <- as.character(clusters$Protein[clusters$Clust2 == cluster])
  proteins <- proteins[!is.na(proteins)]
  probes <- dic_transcriptogramer$PROBE[dic_transcriptogramer$ENSP %in% proteins]
  
  #cores
  cols <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  #subset da matrix de coexpressão contendo proteinas do cluster
  res2 <- ceaCpm[rownames(ceaCpm) %in% probes, colnames(ceaCpm) %in% probes ]
  #grapho a partir da matriz de adjacência
  ceg <- graph.adjacency(res2, diag=FALSE, mode="undirected", weighted=TRUE)
  #gt<-ceg
  gt <- subg(g=ceg, dat=DElimma2, refcol=8)
  gt <- att.setv(g=gt, from="Symbol", to="nodeAlias")
  #gt <- att.setv(g=gt, from=logfc_grupo, to="nodeColor", breaks = seq(-5,5,0.5), pal=2)
  gt <- att.setv(g=gt, from="logFC", to="nodeColor", cols = cols, breaks = breaks)
  resetd(rdp)
  E(gt)$edgeWidth <- 1
  E(gt)$edgeColor <- "#d1d1d1ff"
  V(gt)$nodeShape <- ifelse(V(gt)$nodeAlias %in% tfs, "DIAMOND", "CIRCLE")
  V(gt)$nodeSize <- ifelse(V(gt)$nodeAlias %in% tfs, 60, 30)
  
  V(gt)$nodeLineColor <- "#dadadaff"
  V(gt)$nodeLineWidth <- 0
  V(gt)$nodeAlias <- ""
  addGraph(rdp, gt, gcoord=c(10,25), gscale=50, theme='tm1', zoom=100)
  hc <- hclust(dist(get.adjacency(gt, attr="weight",sparse=F)))
  nesthc(rdp,hc, cutlevel=cut_level, nmemb=5, cex=0.3, labels=V(gt)$nodeAlias)
  mergeOutEdges(rdp,nlev=2)
  relax(rdp,25,10,50,150,100,75,5,20,ps=F)
  
  scl <- gt$legNodeColor$scale
  leg <- gt$legNodeColor$legend
  addLegend.color(rdp, colvec=scl, labvec=leg, title = "Log2 FC")
  
}

cluster=2
logfc_grupo="logFC.Lead30_5...Control_5"
cut_level = 2
breaks = seq(-5, 5, 0.5)

plot_nests_names <- function(cluster, logfc_grupo, cut_level = 2, 
                             breaks = seq(-5, 5, 0.5)) {
  
  print(paste0("Rede - cluster ", cluster))
  
  proteins <- as.character(clusters$Protein[clusters$Clust2 == cluster])
  proteins <- proteins[!is.na(proteins)]
  probes <- dic_transcriptogramer$PROBE[dic_transcriptogramer$ENSP %in% proteins]
  
  
  cols <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  res2 <- ceaCpm[rownames(ceaCpm) %in% probes, colnames(ceaCpm) %in% probes ]
  ceg <- graph.adjacency(res2, diag=FALSE, mode="undirected", weighted=TRUE)
  gt <- subg(g=ceg, dat=DElimma[DElimma$degenes.Lead30_5...Control_5 != 0, ], refcol=3)
  gt <- att.setv(g=gt, from="hgnc_symbol", to="nodeAlias")
  #gt <- att.setv(g=gt, from=logfc_grupo, to="nodeColor", breaks = seq(-5,5,0.5), pal=2)
  gt <- att.setv(g=gt, from=logfc_grupo, to="nodeColor", cols = cols, breaks = breaks)
  resetd(rdp)
  E(gt)$edgeWidth <- 1
  E(gt)$edgeColor <- "#d1d1d1ff"
  V(gt)$nodeShape <- ifelse(V(gt)$nodeAlias %in% tfs, "DIAMOND", "CIRCLE")
  V(gt)$nodeSize <- ifelse(V(gt)$nodeAlias %in% tfs, 60, 30)
  
  V(gt)$nodeLineColor <- "#dadadaff"
  V(gt)$nodeLineWidth <- 0
  addGraph(rdp, gt, gcoord=c(10,25), gscale=50, theme='tm1', zoom=100)
  hc <- hclust(dist(get.adjacency(gt, attr="weight",sparse=F)))
  pdf(file = paste0("Grupo5/clusters", cluster),
      width = 11, height = 7)
  nesthc(rdp,hc, cutlevel=cut_level, nmemb=5, cex=0.3, labels=V(gt)$nodeAlias)
  dev.off()
  mergeOutEdges(rdp,nlev=2)
  relax(rdp,25,10,50,150,100,75,5,20,ps=F)
  
  scl <- gt$legNodeColor$scale
  leg <- gt$legNodeColor$legend
  addLegend.color(rdp, colvec=scl, labvec=leg, title = "Log2 FC")
  
  
}



### LOGFC DO GRUPO 5 --------------------------------------------------------
### redes e dendrogramas serao salvos na pasta 'Grupo5'

plot_nests(cluster = 1, logfc_grupo = "logFC.Lead30_5...Control_5")
plot_nests_names(cluster = 2, logfc_grupo = "logFC.Lead30_5...Control_5")

plot_nests(cluster = 3, logfc_grupo = "logFC.Lead30_5...Control_5")
plot_nests_names(cluster = 3, logfc_grupo = "logFC.Lead30_5...Control_5")

plot_nests(cluster = 6, logfc_grupo = "logFC.Lead30_5...Control_5")
plot_nests_names(cluster = 5, logfc_grupo = "logFC.Lead30_5...Control_5")

plot_nests(cluster = 7, logfc_grupo = "logFC.Lead30_5...Control_5")
plot_nests_names(cluster = 6, logfc_grupo = "logFC.Lead30_5...Control_5")

plot_nests(cluster = 8, logfc_grupo = "logFC.Lead30_5...Control_5")
plot_nests_names(cluster = 8, logfc_grupo = "logFC.Lead30_5...Control_5")

plot_nests(cluster = 9, logfc_grupo = "logFC.Lead30_5...Control_5")
plot_nests_names(cluster = 9, logfc_grupo = "logFC.Lead30_5...Control_5")

plot_nests(cluster = 10, logfc_grupo = "logFC.Lead30_5...Control_5")
plot_nests_names(cluster = 16, logfc_grupo = "logFC.Lead30_5...Control_5")

plot_nests(cluster = 12, logfc_grupo = "logFC.Lead30_5...Control_5")
plot_nests_names(cluster = 17, logfc_grupo = "logFC.Lead30_5...Control_5")





