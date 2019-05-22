rm(list = ls())

library(affy)
library(factoextra)

setwd("/home/clovis/Dropbox/Chumbo/PCA/")

full_datasets_path <- "./datasets_chumbo/"
dts <- dir(full_datasets_path)

dataset="GSE28261"
norm_and_pca <- function(dataset) {
  # pch.code<-c(20,4)
  raw <- ReadAffy(celfile.path = paste0(full_datasets_path, dataset))
  norm <- rma(raw)
  my_exp <- exprs(norm)
  
  samples <- read.table(paste0(full_datasets_path, dataset, "/samples.txt"), 
                        header = F, stringsAsFactors = F, sep = "\t", 
                        col.names = c("GSMid", "groups", "type"))
  
  targets <- as.factor(paste(samples$groups,samples$type))
  data.exprs <- as.data.frame(na.omit(my_exp))
  colnames(data.exprs) <- targets
  tdata.exprs <-t(data.exprs)
#  color.code <- colorRampPalette(c("blue", "red", "green", "orange", "brown"))(length(levels(targets)))
  pca <- prcomp(tdata.exprs)
  #calcula os eingvectors e proporção explicada
  tmp <- pca$sdev^2
  eigs<-tmp / sum(tmp)
  
  pcs<-data.frame(pca$x[,1:2])
  # pcs$color<-"0"
  pcs$shape<-"0"
  names<-rownames(pcs)
  names<-gsub(pattern = "[.]ctrl.*",replacement = "", names)
  names<-gsub(pattern = "[.]exp.*",replacement = "", names)
  names<-gsub(pattern = "[.]",replacement = "_", names)
  names<-gsub(pattern = "_1",replacement = "", names)
  names<-gsub(pattern = "_2",replacement = "", names)
  pcs$name<-names
  controls<-grep(pattern = "ctrl", rownames(pcs))
  exposure<-grep(pattern = "exp", rownames(pcs))
  pcs$shape[controls]<-"20"
  pcs$shape[exposure]<- "4"
  
  names<-data.frame(unique(names))
  names$colors<-colorRampPalette(c("blue","darkgreen","red","deeppink", "darkorchid", "orange"))(nrow(names))
  colnames(names)<-c("name", "colors")
  
  pcs<-merge(pcs, names, by="name")
  
#  expPalet<-colorRampPalette(c("red","deeppink", "gold","darkorchid"))(length(exposure))
#  pcs$color[exposure]
  
  pcCtrl<-pcs[pcs$shape == 20,]
  pcCase<-pcs[pcs$shape == 4,]
  ctColors<-pcCtrl$color
  csColors<-pcCase$color
  ctNames<-pcCtrl$name
  csNames<-pcCase$name
  
 # colorName<-unique(pcs[,c(1,5)])
  
  g<-ggplot()+
    ggtitle(dataset)+
    xlab(paste0("PC1 (",round(eigs[1]*100,digits = 1),"%)"))+
    ylab(paste0("PC2 (",round(eigs[2]*100,digits = 1),"%)"))+
    geom_point(data = pcCtrl,
               aes(PC1,PC2,
                   #label=names,
                   shape=rep("20",length(ctColors)),
                   color = ctColors))+
    geom_point(data = pcCase,
               aes(PC1,PC2,
                   #label=names,
                   shape=rep("4",length(csColors)),
                   color = csColors))+
    theme_bw()+
    scale_shape_manual(name = "",
                       guide = "legend",
                       values = c(20,4),
                       labels=c("Control","Treated"))+
    scale_color_manual(name = "",
                       guide = "legend",
                       values = names$colors,
                       labels=names$name)+
    
    theme(legend.position = "right",
          legend.background = element_rect(fill = alpha("white",0)))
  
  pdf(file = paste0("/home/clovis/Dropbox/Chumbo/PCA/PCA_", dataset,".pdf"), 
      width = 8, height = 5.8)
  plot(g)
  dev.off()
  
  # pdf(file = paste0(full_datasets_path, dataset, "/PCA_", dataset), 
  #     width = 10, height = 7)
  # plot(x = pca$x[,1:2], pch = 19, col = color.code[targets], main = "PCA after normalization")
  # text(pca$x[,1]-1, pca$x[,2]-1, samples$groups, cex = 0.5) 
  # dev.off()
  
  
}



coefVar = function(samples, dataset){
  raw <- ReadAffy(celfile.path = paste0(full_datasets_path, dataset))
  norm <- rma(raw)
  my_exp <- exprs(norm)
  samples <- read.table(paste0(full_datasets_path, dataset, "/samples.txt"), 
                        header = F, stringsAsFactors = F, sep = "\t", 
                        col.names = c("GSMid", "groups"))
  
  sNames<-unique(samples[,2])
  name="C57BL_6J_exp"
  log<-paste0(full_datasets_path, dataset, "/result.txt")
  for(name in sNames){
    columns<-paste0(samples$GSMid[samples$groups == name],".CEL.gz")
    cat("Sample ",name,": ",length(columns),"replic.\n",
        file = log,append = T )
    # if(length(columns)==1){
    #   next
    # }
    tmp<-data.frame(my_exp[,colnames(my_exp)%in%columns])
    tmp$mean<-(apply(X = tmp,MARGIN = 1,FUN = mean))
    tmp$sd<-apply(X = tmp,MARGIN = 1,FUN = sd)
    tmp$CV<-tmp$sd/tmp$mean
    cat("\tMean CV:",mean(tmp$CV),".\n",
        file = log,append = T)
  }
}

lapply(dts, norm_and_pca)

