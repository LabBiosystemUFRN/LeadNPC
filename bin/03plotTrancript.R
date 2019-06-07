rm(list = ls())

figuras="figuras"

setwd("/home/clovis/Dropbox/Chumbo/")
#load("./Data/pheno_data.RData")
load("./Data/counts.RData")
load("./transcriptograms/allTranscriptogramers80")

#check de transc
transc[[1]]@clusters
transc[[2]]@clusters

###########################
#cores sem bound
rainbowcols<-c("#FF2222FF","#ff9e30ff","#00AAAA77","#ffff00ff","#8000AA77",
               "#FF00BF77","#6666ff77","#80660066","#FF222277","#00FF40FF",
               "#ff9e3077","#ffff0077","#00FF4077","#00AAAAFF","#0000ffff",
               "#00550077","#00005577","#000055ff","#AA0000FF","#bf6e00ff",
               "#005500FF","#80660066","#999900ff","#00FFFFFF","#8000AAFF",
               "#3333aaff","#6666ffff","#FFffBFFF","#FF00BFFF")
# barplot(c(seq(1:29)),
#         col = rainbowcols,names.arg = seq(1:29))
# rainbowcols2<-c("#FF222277","#ff9e3077","#ffff0077",
#                 "#00FF4077","#0000ff77","#80660066",
#                 "#AA000077","#bf6e0077",,
#                 "#99990077","#00FFFF77","#00AAAA77","#8000AA77",
#                 "#FF00BF77","#6666ff77","#3333aa77")
# barplot(c(seq(1:18)),
#         col = rainbowcols2,names.arg = seq(1:18))

color<-list()
color[[1]]<-rainbowcols[c(1,2,4,10,11,14,15,18,19,20,21,22,23,25,29)]
color[[2]]<-rainbowcols[c(1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,24,25,26,27,28,29)]
#save(color,file = "./Data/colors.RData")

#"#3333aa77"

#############################
#cores anteriores
# rainbowcols<-c("#FF2222FF","#ff9e30ff","#ffff00ff","#009900FF",
#                "#00FF40FF","#0000ffff","#80660066",
#                "#AA0000FF","#bf6e00ff","#005500FF","#000055ff",
#                "#999900ff","#00FFFFFF","#00AAAAFF","#8000AAFF",
#                "#FF00BFFF","#6666ffff","#3333aaff")
# rainbowcols2<-c("#FF222277","#ff9e3077","#ffff0077","#00990077",
#                "#00FF4077","#0000ff77","#80660066",
#                "#AA000077","#bf6e0077","#00550077","#00005577",
#                "#99990077","#00FFFF77","#00AAAA77","#8000AA77",
#                "#FF00BF77","#6666ff77","#3333aa77")
# 
# color<-list()
# 
# 
# 
# 
# color[[1]]<-rainbowcols[c(1	,4	,3	,2  		,6	,7	,11	,14	,13	,15	,16   			,1)]
# color[[2]]<-rainbowcols[c(1   			,2	,5	,6	,7	,11	,14	,15 		,16	,17	,18	,1)]
#save(color,file = "./Data/colors.RData")
# color[[3]]<-rainbowcols[c(2,4,9,11,14,17,22,18)]
# color[[4]]<-rainbowcols[c(1,5,2,3,4,6,10,8,9,11,7,14,17,12,15,16,13,22,18,23,21,20,19,24,25)]
# color[[5]]<-rainbowcols[c(1,2,4,8,9,11,7,14,6,10,12,21,15,16,13,22,18,23,20,24,21,19,23,25)]
# color[[1]]<-rainbowcols[c(1,4,9,10,14,19,22,13,2,)]

# rainbowcols <- rainbow(8)
# barplot(c(seq(1:18)),
#         col = rainbowcols,names.arg = seq(1,18))
# barplot(c(seq(1:11)),
#         col = color[[2]],names.arg = c(1,2,3,5,6,7,11,14,15,16,17))
# barplot(c(seq(2,25,3)),
#         col = rainbowcols[seq(2,25,3)],names.arg = seq(1,25,3))

lim<-0
i=2
setwd("./levels/")
for(i in 1:2){
  object <- transc[[i]]
  
  
  levels<-read.table(file = paste0("AllW80LevelsPCAGr",
                                   i,"Lead30.csv")
                     ,sep ="\t",stringsAsFactors = F)
  levels<-as.logical(unname(unlist(levels[1,])))
  
  
  case <- object@transcriptogramS2[, -c(1, 2)]
  control <- case[, which(levels == TRUE)]
  case <- case[, which(levels == FALSE)]
  n <- nrow(control)
  caseValues <- vapply(seq.int(1, n), function(i) {
    result <- mean(unlist(case[i, ])) - mean(unlist(control[i,
                                                            ]))
    return(result)
  }, numeric(1))
  
  
  cat(max(caseValues),"\n")
  pBreaks<- object@clusters
  smoothedLine <- stats::smooth.spline(object@transcriptogramS2$Position,
                                       caseValues, spar = 0.35)
  #lim <- max(abs(caseValues))
  rm(case, control, n, caseValues)
  #lim <- round(lim, digits = 1)
  lim<-0.65
  cat(lim)
  if(object@pbc){
    myColors <- color[[i]]
    myColors <- c(myColors, myColors[1])
  }else{
    myColors <- color[[i]]
  }
  df <- data.frame(x = smoothedLine$x, y = smoothedLine$y)
  rm(smoothedLine)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_line(lwd = 1, ggplot2::aes(y = 0, colour = "c1")) +
    ggplot2::geom_line(lwd = 1, ggplot2::aes(y = y, colour = "c2")) +
    ggplot2::scale_y_continuous(limits = c(-lim, lim), breaks = seq(-lim, lim, 0.1)) +
    ggplot2::scale_x_continuous(limits = c(0, length(object@ordering$Position) - 1),
                                breaks = seq.int(0, length(object@ordering$Position) - 1, 1000)) +
    ggplot2::scale_colour_manual(values = c("black", "grey80"), name = "Conditions",
                                 labels =  c("Control", "Case")) +
    ggplot2::scale_linetype_manual(values = "blank", name = "Number of clusters",
                                   labels = ifelse(object@pbc, length(myColors) - 1, length(myColors))) +
    ggplot2::labs(x = "Gene position", y = "Difference of means (case - control)") +
    ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  invisible(sapply(seq.int(1, length(pBreaks)), function(i) {
    idx <- which(df$x >= pBreaks[[i]][1] & df$x <= pBreaks[[i]][2])
    aux <- data.frame(x = df$x[idx], y = df$y[idx])
    p <<- p + ggplot2::geom_line(data = aux, lwd = 1.4, col = myColors[i],
                                 ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line(ggplot2::aes(linetype = "lines"))
    return(NULL)
  }))

  plot(p)  
  pdf(width = 11,height = 4.75,file = paste0("../",figuras,"/transcGrupo",i,".pdf"))
  suppressMessages(graphics::plot(p))
  dev.off()
}



# color<-list()
# color[[1]]<-rainbowcols2[c(2,3,4,5,6,8,14,15,16,9,10,11,12,13,17,18,1)]
# color[[2]]<-rainbowcols2[c(1,2,3,5,6,8,7,11,14,15,16,17)]

lim<-0
i=2
#setwd("/home/clovis/Dropbox/Clovis/analise/Novo/levels/")
p<-ggplot2::ggplot()
for(i in 1:2){
  object <- transc[[i]]
  
  
  levels<-read.table(file = paste0("AllW125LevelsPCAGr",
                                   i,"Lead30.csv")
                     ,sep ="\t",stringsAsFactors = F)
  levels<-as.logical(unname(unlist(levels[1,])))
  
  
  case <- object@transcriptogramS2[, -c(1, 2)]
  control <- case[, which(levels == TRUE)]
  case <- case[, which(levels == FALSE)]
  n <- nrow(control)
  caseValues <- vapply(seq.int(1, n), function(i) {
    result <- mean(unlist(case[i, ])) - mean(unlist(control[i,
                                                            ]))
    return(result)
  }, numeric(1))
  
  
  cat(max(caseValues),"\n")
  pBreaks<- object@clusters
  smoothedLine <- stats::smooth.spline(object@transcriptogramS2$Position,
                                       caseValues, spar = 0.35)
  #lim <- max(abs(caseValues))
  rm(case, control, n, caseValues)
  #lim <- round(lim, digits = 1)
  lim<-0.4
  cat(lim)
  if(object@pbc){
    myColors <- color[[i]]
    myColors <- c(myColors, myColors[1])
  }else{
    myColors <- color[[i]]
  }
  df <- data.frame(x = smoothedLine$x, y = smoothedLine$y)
  rm(smoothedLine)
 # p<-ggplot2::ggplot()
  p <-p  +
    ggplot2::geom_line(data = df,lwd = 1, ggplot2::aes(x=x,y = 0, colour = "c1")) +
    ggplot2::geom_line(data = df,lwd = 1, ggplot2::aes(x=x,y = y, colour = "c2")) +
    ggplot2::scale_y_continuous(limits = c(-lim, lim), breaks = round(seq(-lim, lim, 0.1), digits = 1)) +
    ggplot2::scale_x_continuous(limits = c(0, length(object@ordering$Position) - 1),
                                breaks = seq.int(0, length(object@ordering$Position) - 1, 1000)) +
    ggplot2::scale_colour_manual(values = c("black", "grey80"), name = "Conditions",
                                 labels =  c("Control", "Case")) +
    ggplot2::scale_linetype_manual(values = "blank", name = "Number of clusters",
                                   labels = ifelse(object@pbc, length(myColors) - 1, length(myColors))) +
    ggplot2::labs(x = "Gene position", y = "Difference of means (case - control)") +
    ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  invisible(sapply(seq.int(1, length(pBreaks)), function(i) {
    idx <- which(df$x >= pBreaks[[i]][1] & df$x <= pBreaks[[i]][2])
    aux <- data.frame(x = df$x[idx], y = df$y[idx])
    p <<- p + ggplot2::geom_line(data = aux, lwd = 1.4, col = myColors[i],
                                 ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line(data = df,ggplot2::aes(x=x,y=y,linetype = "lines"))
    return(NULL)
  }))
  
}
pdf(width = 11,height = 7,file = paste0("../",figuras,"/transcAmbos.pdf"))
suppressMessages(graphics::plot(p))
dev.off()



smoothedLine<-list()
lim<-0
for(i in 1:2){
  object <- transc[[i]]
  
  
  levels<-read.table(file = paste0("AllW125LevelsPCAGr",
                                   i,"Lead30.csv")
                     ,sep ="\t",stringsAsFactors = F)
  levels<-as.logical(unname(unlist(levels[1,])))
  
  
  case <- object@transcriptogramS2[, -c(1, 2)]
  control <- case[, which(levels == TRUE)]
  case <- case[, which(levels == FALSE)]
  n <- nrow(control)
  caseValues <- vapply(seq.int(1, n), function(i) {
    result <- mean(unlist(case[i, ])) - mean(unlist(control[i,
                                                            ]))
    return(result)
  }, numeric(1))
  
  
  cat(max(caseValues),"\n")
  pBreaks<- object@clusters
  smoothedLine[[i]] <- stats::smooth.spline(object@transcriptogramS2$Position,
                                       caseValues, spar = 0.35)
}


  if(object@pbc){
    myColors <- color[[2]]
    myColors <- c(myColors, myColors[2])
  }else{
    myColors <- color[[2]]
  }
  df <- data.frame(x = smoothedLine[[2]]$x, y = (smoothedLine[[2]]$y-smoothedLine[[1]]$y))
  rm(smoothedLine)
  #lim <- max(abs(df$y))
 # lim <- round(lim, digits = 1)
  lim<-0.4
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_line(lwd = 1, ggplot2::aes(y = 0, colour = "c1")) +
    ggplot2::geom_line(lwd = 1, ggplot2::aes(y = y, colour = "c2")) +
    ggplot2::scale_y_continuous(limits = c(-lim, lim), breaks = round(seq(-lim, lim, 0.1), digits = 1)) +
    ggplot2::scale_x_continuous(limits = c(0, length(object@ordering$Position) - 1),
                                breaks = seq.int(0, length(object@ordering$Position) - 1, 1000)) +
    ggplot2::scale_colour_manual(values = c("black", "grey80"), name = "Conditions",
                                 labels =  c("Control", "Case")) +
    ggplot2::scale_linetype_manual(values = "blank", name = "Number of clusters",
                                   labels = ifelse(object@pbc, length(myColors) - 1, length(myColors))) +
    ggplot2::labs(x = "Gene position", y = "Difference of means (case - control)") +
    ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  invisible(sapply(seq.int(1, length(pBreaks)), function(i) {
    idx <- which(df$x >= pBreaks[[i]][1] & df$x <= pBreaks[[i]][2])
    aux <- data.frame(x = df$x[idx], y = df$y[idx])
    p <<- p + ggplot2::geom_line(data = aux, lwd = 1.4, col = myColors[i],
                                 ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line(ggplot2::aes(linetype = "lines"))
    return(NULL)
  }))
  
  pdf(width = 11,height = 7,file = paste0("../",figuras,"/transcDif.pdf"))
  suppressMessages(graphics::plot(p))
  dev.off()






# retorno<-data.frame(x=numeric(),y=numeric(),cor=numeric(),stringsAsFactors = F)
# i=2
# for( i in seq.int(1, length(pBreaks))){
#   idx <- which(df$x >= pBreaks[[i]][1] & df$x <= pBreaks[[i]][2])
#   aux <- data.frame(x = df$x[idx], 
#                     y = df$y[idx],
#                     cor=rep(color[[2]][i],length(idx)),
#                     stringsAsFactors = F)
#   retorno<-rbind(retorno,aux)
# }
# retorno$cor[is.na(retorno$cor)]<-"#000000FF"
# ggplot2::ggplot()+
#   ggplot2::geom_line(data = retorno, lwd = 1.4, 
#                      ggplot2::aes(x = x, y = y, color = cor))

  
