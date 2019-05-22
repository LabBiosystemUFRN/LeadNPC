rm(list = ls())

library(ggplot2)
library(biomaRt)
library(hgu95av2.db)
library(topGO)

#calcula a primeira derivada do gráfico
fderivada = function(df){
  derivada <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(derivada)<-c("x","y")
  for(i in seq(1,length(df$x)-1,1)){
    #janela <- df[(i-2):(i+2),]
    derivada[i,2] <- (df[i,2]-df[i+1,2])/
      (df[i,1]-df[i+1,1])
    derivada[i,1] <- df[i,1]
  }
  return(derivada)
}


#localiza picos do gráfico
achaPico = function (derivada){
  #compara os sinais da derivada. Havendo diferença houve inflexão no gráfico
  sinal<- sign(derivada$y)
  #variavel de resultado
  resultado<-data.frame(matrix(ncol = 2, nrow = 0))
  tmp<-data.frame(matrix(ncol = 2, nrow = 1))
  for(i in seq(1,length(sinal)-1,1)){
    if(sinal[i] != sinal[i+1]){
      #em tipo picos recebem 1, vales recebem -1
      if(sinal[i] > sinal[i+1]){
        tmp[1,2]<-1
      }else{
        tmp[1,2]<- -1
      }
      tmp[1,1]<-derivada$x[i]
      
      resultado<- rbind(resultado,tmp)
    }
  }
  colnames(resultado)<-c("x","tipo")
  return(resultado)
}

#corta picos e vales que não ultrapassem o threshold
aplicaTH= function (df, picos, limiarUp){
  #Cancelar picos e vales menores q limiarDw
  #Deve haver N picos e N+1 vales. Acrescentar no inicio e fim se preciso
  # Como valores são normalizados entre 0 e 1 colocar primeiro vale em -0.1 e ultimo em 1
  #Deve haver alternancia entre picos e vales tb. caso não haja, acusar erro de processamento
  #testar se pico é excede os vales adjacntes em limiarUp
  #caso contrário, cortar pico e vale posterior
  
  #cria variavel de retorno
  resultado<-data.frame(matrix(ncol = 3, nrow = 0))
  #df temporario juntando tabela original e picos
  tmp <- merge(df,picos, by= 1)
  colnames(tmp)<- c("x","y","tipo")
  #erro se tmp ficar vazio
  if(nrow(tmp)<=1){
    cat(paste("Nenhum pico e vale a ser processado... \n"))
    #print("Erro: ")
    colnames(resultado)<-c("x","y","tipo")
    return(resultado)
  }
  
  #acrescenta o primeiro vale, caso este não exista
  if(tmp$tipo[1] == 1){
    tmp<-rbind(c(-0.1,-1,-1),tmp)
  }
  #acrescenta o ultimo vale, caso este não exista
  if(tmp$tipo[nrow(tmp)] == 1){
    tmp<-rbind(tmp,c(max((tmp$x)+1000),-1,-1))
  }
  #verifica se existe alternancia entre picos e vales
  for(i in seq(1,nrow(tmp)-1,2)){
    #primeiro vale depois pico
    if(tmp$tipo[i]== 1 | tmp$tipo[i+1] == -1 ){
      cat(paste("Inconsistência no número de picos e vales... \n"))
      #      print("Erro: ")
      colnames(resultado)<-c("x","y","tipo")
      return(resultado)
    }
  }

  while(nrow(tmp) >= 3){
    #realiza verificação onde vale/pico/vales > limiarUp
    valeAnt<-tmp$y[1]+1
    picoAnt<-tmp$y[2]+1
    valePos<-tmp$y[3]+1
    if(nrow(tmp) >= 4){
      picoPos<-tmp$y[4]+1
    }else{
      picoPos<- 1
    }
    if(limiarUp < picoAnt - valeAnt & limiarUp < picoAnt - valePos){
      resultado <- rbind(resultado, tmp[2,],tmp[3,])
    }
    if(limiarUp < picoAnt - valePos & limiarUp < picoPos - valePos){
      resultado <- rbind(resultado, tmp[3,],tmp[4,])
    }
    tmp <- tmp[-c(1,2),]
    # }else{
    #   tmp <- tmp[-c(2,3),]
    # }
  }
  colnames(resultado)<-c("x","y","tipo")
  resultado<-unique(resultado)
  return(resultado)
}

sigGOs = function(myInterestedGenes,allGenes,pval){
  geneList <- factor(as.integer(allGenes %in% myInterestedGenes))
  names(geneList) <- allGenes
  
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                nodeSize = 5,
                annot = annFUN.org,
                mapping = "org.Hs.eg.db",
                ID = "ensembl")
  
  sg <- sigGenes(GOdata)
  numSigGenes(GOdata)
  resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher") 
  
  allRes <- GenTable(GOdata, 
                     classicFisher = resultFisher, 
                     orderBy = "resultFisher", 
                     ranksOf = "classicFisher",
                     topNodes = length(usedGO(GOdata)))
  
  allRes$pAdj<-p.adjust(p = allRes$classicFisher,method = "BH")
  sigGo<-allRes[allRes$pAdj<=pval,]
  return(sigGo)
}

# Principal ----
figuras="figuras"

setwd("/home/clovis/Dropbox/Chumbo/")
load("./Data/counts.RData")
load("./transcriptograms/allTranscriptogramers80")
load(file = "./Data/colors.RData")
load(file= "./Data/dictPep2Gene.RData")

#gera lista de ensomblId por posição no transcriptograma
genePos<-transc[[1]]@ordering
genePos<-merge(genePos,dic_transcriptogramer,
               by.x="Protein", 
               by.y="ENSP")
colnames(genePos)<-c("ProteinID","Position","EnsemblID")

#gera lista de todos os genes
xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "ensembl")
head(xx)
allGenes <- unique(unlist(xx))

radius<-80
lim<-0
i=1
setwd("./levels/")
limiarUp = 0.025
for(i in 1:1){
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
  lim<-0.65+0.025
  cat(lim)
  if(object@pbc){
    myColors <- color[[i]]
    myColors <- c(myColors, myColors[1])
  }else{
    myColors <- color[[i]]
  }
  df <- data.frame(x = smoothedLine$x, y = smoothedLine$y)
  rm(smoothedLine)
  
  derivada<-fderivada(df)
  picos<-achaPico(derivada)
  picos<-aplicaTH(df = df,picos = picos,limiarUp = limiarUp )
  tmp<-picos[1,]
  tmp$cluster<-0
  tmp<-tmp[-1,]
  cl=3
  #informações sobre os limites dos clusteres
  limClst<-data.frame(cluster=integer(),
                      inicio=integer(),
                      fim = integer())
  for(cl in 1:length(transc[[i]]@clusters)){
    inicio=transc[[i]]@clusters[[cl]][1]
    fim=transc[[i]]@clusters[[cl]][2]
    tmp2<-picos[picos$x>=inicio & picos$x<=fim,]
    tmp2$cluster<-cl
    
    tmp<-rbind(tmp,tmp2 )
    limClst[nrow(limClst)+1,]<-c(cl,inicio,fim)
  }
  
  picos<-na.exclude(tmp)
  picos$nr<-seq(1,nrow(picos))
  
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
  sigGo<-integer()
  nr=1
  for(nr in 1:nrow(picos)){
    p<-p+annotate(geom="text", 
                  x=picos$x[nr], 
                  y=picos$y[nr] +(0.025*picos$tipo[nr]), 
           label=picos$nr[nr],
           color="black",cex=2)
    #pega os proteinId dentro da janela
    #testa se janela ultrapassa limites 
    if(picos$x[nr]>=radius){
      inicio=(picos$x[nr]-radius)
    }else{
      inicio=0
    }
    if(max(genePos$Position)-picos$x[nr]>=radius){
      fim=(picos$x[nr]+radius)
    }else{
      fim=max(genePos$Position)
    }
    #remove genes fora do cluster
    if(inicio<limClst$inicio[picos$cluster[nr]]){
      inicio<-limClst$inicio[picos$cluster[nr]]
    }
    if(fim>limClst$fim[picos$cluster[nr]]){
      fim<-limClst$fim[picos$cluster[nr]]
    }
    myGenes<-unique(genePos$EnsemblID[genePos$Position>=inicio&
                         genePos$Position<=fim])
    cat("Enriquecimento tempo", i, "pico", nr,"de",nrow(picos),"\n")
    tmp<- sigGOs(myInterestedGenes = myGenes,
                  allGenes = allGenes, 
                  pval = 0.01)
    tmp$Cluster<-picos$cluster[nr]
    tmp$peak<-nr
    tmp<-tmp[,c("Cluster","peak","GO.ID","Term","Annotated","Significant","Expected",   
                    "classicFisher","pAdj")]
    if(typeof(sigGo)=="integer"){
      sigGo<-tmp
    }else{
      sigGo<-rbind(sigGo,tmp)
    }
  }
  # plot(p) 
  write.csv(sigGo,file = paste0("../terms/PicosW80PCAGr",i,".csv"),row.names = F)
  pdf(width = 11,height = 4.75,file = paste0("../",figuras,"/transcGrupo",i,"Picos.pdf"))
  suppressMessages(graphics::plot(p))
  dev.off()
}

