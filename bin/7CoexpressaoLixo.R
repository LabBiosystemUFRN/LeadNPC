
rm(list = ls())

#library(edgeR)
#library(refGenome)
library(transcriptogramer)
library(biomaRt)
library("FactoMineR")
library("factoextra")
library(ggplot2)
library(RedeR)
library(igraph)
library(dplyr)



# Ler tabela de counts ----------------------------------------------------
setwd("/home/clovis/Dropbox/Chumbo/Data/")
load("./allTranscriptogramers80")
load(file = "./counts.RData")
load(file = "./associationHs700.RData")
load(file = "./CPM.RData")

# #Verifica se todas as mostras relevantes estão presentes
# n01<-data.frame(run=pheno_data$Run[pheno_data$source_name == "Lead30_NPCs"],
#                 treat=pheno_data$source_name[pheno_data$source_name == "Lead30_NPCs"],
#                 day=pheno_data$Day[pheno_data$source_name == "Lead30_NPCs"],
#                 stringsAsFactors = F)
# n02<-data.frame(run=pheno_data$Run[pheno_data$source_name == "Control_NPCs"],
#                 treat=pheno_data$source_name[pheno_data$source_name == "Control_NPCs"],
#                 day=pheno_data$Day[pheno_data$source_name == "Control_NPCs"],
#                 stringsAsFactors = F)
# #n2<-data.frame(run=colnames(DEexp),stringsAsFactors = F)
# #n2$id<-"n2"
# n3<-rbind(n01,n02)
# #todas<- merge(n2,n3,by="run",all = T)
# #rm(list = c("n01","n02","n1","n2","n3","todas"))

# var set ----
name = "Lead30"
set<-"All"
radius<-0
pval <- 0.001
transc<-list()

#calculando expressão diferencial de todas as amostras
group<-c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
case_control30 <- pheno_data$Run[pheno_data$source_name == "Control_NPCs"&
                                     pheno_data$Day%in%group]
  
  
  
case_control30 <- append(case_control30, pheno_data$Run[pheno_data$source_name == paste0(name,"_NPCs")&
                                                            pheno_data$Day%in%group])
expression<- as.data.frame(logCPM[, colnames(logCPM)%in%case_control30])

t1 <- transcriptogramPreprocess(association = association, 
                                  ordering = Hs700,
                                  radius = radius )
  
t1 <- transcriptogramStep1(object = t1, expression = expression,
                             dictionary = dic_transcriptogramer, nCores = T)
  
t1 <- transcriptogramStep2(object = t1, nCores = T)
  
levels <- c(rep(TRUE,length(group)),rep(FALSE,length(group)))
tmp2<-rbind(levels,case_control30)

t1 <- differentiallyExpressed2(object = t1, 
                              levels = levels, 
                              pValue = pval, 
                              species = "Homo sapiens",
                              boundaryConditions = T)

#save(t1,file = "t1.RData")


check_species1 <- function(argument){
  if(is.data.frame(argument)){
    colnames(argument) <- c("ensembl_peptide_id", "external_gene_name")
    if(species1Check(argument)){
      return (argument)
    }else{
      stop("argument species - does not have a valid value!",
           call. = FALSE)
    }
  } else if (!((is.character(argument) &&
                length(strsplit(argument, " ")[[1]]) ==
                2 && length(argument) ==
                1))) {
    stop("argument species - does not have a valid value!",
         call. = FALSE)
  }
}

#DE ----
differentiallyExpressed2 = function(object,
                                   levels, pValue = 0.05, species = object@Protein2Symbol, adjustMethod = "BH",
                                   trend = FALSE, title = "Differential expression",
                                   boundaryConditions = FALSE) {
  if (object@status < 2L) {
    stop("argument of class Transcriptogram - be sure to ",
         "call the methods transcriptogramStep1() and ",
         "transcriptogramStep2() before this one!")
  }
  aux <- species
  if (is.data.frame(aux)) {
    species <- check_species1(species)
  } else {
    check_species1(species)
  }
  rm(aux)
  # check_adjustMethod1(adjustMethod)
  # check_trend(trend)
  # check_boundaryConditions(boundaryConditions)
  # check_levels(levels)
  if (length(levels) != (ncol(object@transcriptogramS2) -
                         2)) {
    stop("argument levels - does not have a valid length!")
  }
  object@pbc = FALSE
  levels <- as.factor(levels)
  design <- stats::model.matrix(~0 + levels)
  contrasts <- "levelsFALSE-levelsTRUE"
  fit <- limma::lmFit(as.matrix(object@transcriptogramS2[,
                                                         -c(1, 2)]), design)
  fit$Protein <- object@transcriptogramS2[, 1]
  fit$Position <- object@transcriptogramS2[, 2]
  contrasts <- limma::makeContrasts(contrasts = contrasts,
                                    levels = design)
  rm(design)
  message("calculating statistics... step 1 of 4")
  ct.fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts), trend = trend)
  res.fit <- limma::decideTests(ct.fit, method = "global",
                                adjust.method = adjustMethod, p.value = pValue)
  temp <- data.frame(Protein = ct.fit$Protein, Position = ct.fit$Position,
                     logFC = ct.fit$coef, pValue = ct.fit$p.value,
                     degenes = as.integer(unclass(res.fit)), stringsAsFactors = FALSE)
  rm(contrasts)
  features <- rowSums(res.fit != 0) > 0
  DElimma <- temp[features, ]
  if (nrow(DElimma) == 0) {
    stop("no differentially expressed protein, ",
         "meeting the p-value requirement, was detected!")
  }
  rm(temp)
  colnames(DElimma)[c(3, 4, 5)] <- c("logFC", "pValue", "DEgenes")
  rownames(DElimma) <- NULL
  message("identifying clusters... step 2 of 4")
  pBreaks <- list()
  positions <- DElimma$Position
  clusterStartIndex <- clusterNumber <- 1
  nextIndex <- NULL
  invisible(sapply(seq.int(1, (length(positions) - 1)), function(i) {
    nextIndex <<- i + 1
    if ((positions[nextIndex] - positions[i]) > object@radius) {
      pBreaks[[clusterNumber]] <<- c(positions[clusterStartIndex],
                                     positions[i])
      clusterStartIndex <<- nextIndex
      clusterNumber <<- clusterNumber + 1
    }
    return(NULL)
  }))
  pBreaks[[clusterNumber]] <- c(positions[clusterStartIndex],
                                positions[nextIndex])
  rm(nextIndex, clusterNumber, clusterStartIndex, positions)
  if(boundaryConditions){
    aux <- c()
    min <- object@ordering$Position[1]
    max <- object@ordering$Position[nrow(object@ordering)]
    aux <- invisible(lapply(seq.int(1, length(pBreaks)), function(i) {
      l1 <- pBreaks[[i]][1] - object@radius
      l2 <- pBreaks[[i]][2] + object@radius
      return(c(l1,l2))
    }))
    elim <- list()
    invisible(lapply(seq.int(1, length(aux)), function(i) {
      if(i == length(aux)){
        elim <<- append(elim, list(c(aux[[i]][1], aux[[i]][2])))
      }else if(aux[[i]][2] >= aux[[i + 1]][1]){
        aux[[i + 1]] <<- c(aux[[i]][1], aux[[i + 1]][2])
      }else{
        elim <<- append(elim, list(c(aux[[i]][1], aux[[i]][2])))
      }
      return(NULL)
    }))
    aux <- elim
    if(aux[[1]][1] < min){
      object@pbc = TRUE
      x <- max + 1 + aux[[1]][1] - min
      if(x <= aux[[length(aux)]][2]){
        aux[[1]][1] <- min
        aux[[length(aux)]][2] <- max
      }else{
        aux[[1]][1] <- min
        aux <- append(aux, list(c(x, max)))
      }
    }else if(aux[[length(aux)]][2] > max){
      object@pbc = TRUE
      x <- aux[[length(aux)]][2] %% max + min - 1
      if(x >= aux[[1]][1]){
        aux[[1]][1] <- min
        aux[[length(aux)]][2] <- max
      }else{
        aux[[length(aux)]][2] <- max
        aux <- append(list(c(min, x)), aux)
      }
    }
    pBreaks <- aux
  }
  DElimma$ClusterNumber <- NA
  invisible(sapply(seq.int(1, length(pBreaks)), function(i) {
    if(i == length(pBreaks) && object@pbc){
      DElimma[which(DElimma$Position >= pBreaks[[i]][1] & DElimma$Position <=
                      pBreaks[[i]][2]), "ClusterNumber"] <<- 1
    }else{
      DElimma[which(DElimma$Position >= pBreaks[[i]][1] & DElimma$Position <=
                      pBreaks[[i]][2]), "ClusterNumber"] <<- i
    }
    return(NULL)
  }))
  DElimma <- DElimma[, c(1, 2, 6, 3, 4, 5)]
  # message("generating plot... step 3 of 4")
  # case <- object@transcriptogramS2[, -c(1, 2)]
  # control <- case[, which(levels == TRUE)]
  # case <- case[, which(levels == FALSE)]
  # n <- nrow(control)
  # caseValues <- vapply(seq.int(1, n), function(i) {
  #   result <- mean(unlist(case[i, ])) - mean(unlist(control[i,
  #                                                           ]))
  #   return(result)
  # }, numeric(1))
  # smoothedLine <- stats::smooth.spline(object@transcriptogramS2$Position,
  #                                      caseValues, spar = 0.35)
  # lim <- max(abs(caseValues))
  # rm(case, control, n, caseValues)
  # lim <- round(lim, digits = 1)
  # if(object@pbc){
  #   myColors <- grDevices::rainbow(length(pBreaks) - 1)
  #   myColors <- c(myColors, myColors[1])
  # }else{
  #   myColors <- grDevices::rainbow(length(pBreaks))
  # }
  # df <- data.frame(x = smoothedLine$x, y = smoothedLine$y)
  # rm(smoothedLine)
  # p <- ggplot2::ggplot(df, ggplot2::aes_string("x", "y")) +
  #   ggplot2::geom_line(lwd = 1, ggplot2::aes_string(y = "0", colour = '"c1"')) +
  #   ggplot2::geom_line(lwd = 1, ggplot2::aes_string(y = "y", colour = '"c2"')) +
  #   ggplot2::scale_y_continuous(limits = c(-lim, lim), breaks = round(seq(-lim, lim, 0.1), digits = 1)) +
  #   ggplot2::scale_x_continuous(limits = c(0, length(object@ordering$Position) - 1),
  #                               breaks = seq.int(0, length(object@ordering$Position) - 1, 1000)) +
  #   ggplot2::scale_colour_manual(values = c("black", "grey80"), name = "Conditions",
  #                                labels =  c("Control", "Case")) +
  #   ggplot2::scale_linetype_manual(values = "blank", name = "Number of clusters",
  #                                  labels = ifelse(object@pbc, length(myColors) - 1, length(myColors))) +
  #   ggplot2::labs(x = "Gene position", y = "Difference of means (case - control)", title = title) +
  #   ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  # invisible(sapply(seq.int(1, length(pBreaks)), function(i) {
  #   idx <- which(df$x >= pBreaks[[i]][1] & df$x <= pBreaks[[i]][2])
  #   aux <- data.frame(x = df$x[idx], y = df$y[idx])
  #   p <<- p + ggplot2::geom_line(data = aux, lwd = 1.4, col = myColors[i],
  #                                ggplot2::aes_string(x = "x", y = "y")) +
  #     ggplot2::geom_line(ggplot2::aes_string(linetype = '"lines"'))
  #   return(NULL)
  # }))
  # suppressMessages(graphics::plot(p))
  symbols <- NULL
  message("translating ENSEMBL Peptide ID to SYMBOL... step 3 of 4")
  taxonomyID <- NULL
  if (grepl("\\.", object@ordering[1, 1])) {
    taxonomyID <- strsplit(object@ordering[1, 1], "\\.")[[1]][1]
    taxonomyID <- paste0(taxonomyID, ".")
  }
  if (is.character(species)) {
    message("** this may take some time...")
    species <- tolower(species)
    species <- gsub("^([[:alpha:]]).* ", "\\1", species)
    species <- paste0(species, "_gene_ensembl")
    ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                dataset = species)
    proteins <- NULL
    if (grepl("\\.", object@ordering[1, 1])) {
      proteins <- sapply(strsplit(object@ordering[, 1], "\\."),
                         "[", 2)
    } else {
      proteins <- object@ordering[, 1]
    }
    symbols <- biomaRt::getBM(filters = "ensembl_peptide_id",
                              attributes = c("ensembl_peptide_id", "external_gene_name"),
                              values = proteins, mart = ensembl)
  } else if (is.data.frame(species)) {
    symbols <- species
    if (grepl("\\.", symbols[1, 1])) {
      symbols$ensembl_peptide_id <- sapply(strsplit(symbols[,
                                                            1], "\\."), "[", 2)
    }
  }
  symbols[symbols == ""] <- NA
  symbols <- stats::na.omit(symbols)
  symbols$ensembl_peptide_id <- paste0(taxonomyID,
                                       symbols$ensembl_peptide_id)
  DElimma$Symbol <- NA_character_
  object@Protein2Symbol = symbols
  invisible(sapply(seq.int(1, nrow(DElimma)), function(i) {
    DElimma$Symbol <<- symbols[match(DElimma[, "Protein"],
                                     symbols[,"ensembl_peptide_id"]),
                               "external_gene_name"]
    return(NULL)
  }))
  if(any(is.na(DElimma$Symbol))){
    idx <- which(is.na(DElimma$Symbol))
    DElimma[idx, "Symbol"] <- DElimma[idx, "Protein"]
  }
  object@status = 3L
  object@DE = DElimma
  object@clusters = pBreaks
  message("done!")
  return(object)
}


t1 <- differentiallyExpressed2(object = t1, 
                               levels = levels, 
                               pValue = pval, 
                               species = "Homo sapiens",
                               boundaryConditions = T)




nrow(t1@DE)



#removendo amostras dos dias 0,1,2
exclui<-pheno_data$Run[pheno_data$Day%in%c(0,1,2)]
colunas<-colnames(DEexp)
DEexp<-DEexp[,!colunas%in%exclui]

#limpando o lixo
lista<-ls()
lista<-lista[!lista%in%c("association","DEexp","dic_transcriptogramer")]
rm(list = lista)
rm("lista")

#lendo outros arquivos necessários
load(file = "../transcriptograms/allTranscriptogramers80")
load(file='clusters12.RData')

DEexp<-as.data.frame(DEexp)
DEexp <- tibble::rownames_to_column(DEexp)
DEexp2 <- merge(DEexp,dic_transcriptogramer,
               by.x = "rowname",
               by.y = "PROBE")
DEexp2<-unique(DEexp2)

teste<-DEexp[1:100,2:15]
c1 <- cea2(x= teste, sig=0.001, nper=10000, plot=T, p.adj.method = "BH")
source("../bin/ceaMCore.R")
c2 <- ceaMCore(x= teste, sig=0.001, nper=10000, plot=T, p.adj.method = "BH")


#save(coe_matrix, file = "/home/clovis/Dropbox/Chumbo/Data/coe_matrix.RData")
RedeR::ptcea()
m<-coe_matrix2
m[[1]][1]
nulldist3<-list()
i=1
x<-matrix(nrow = length(m[[1]]$x), ncol = 10000)
x<-data.frame(x)
y<-matrix(nrow = length(m[[1]]$x), ncol = 10000)
y<-data.frame(y)
for (i in 1:10000) {
  x[,i]<- as.vector(m[[i]][1])
  y[,i]<- as.vector(m[[i]][2])
}
