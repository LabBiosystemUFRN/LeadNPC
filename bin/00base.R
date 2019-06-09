rm(list = ls())
baseDir="~/LeadTest/"
setwd(baseDir)

#packages required:
#"biomaRt","affy", "classInt", "data.table","dplyr","edgeR",
#"factoextra","FactoMineR","ggplot2","ggrepel","grid",
#"gridExtra","hgu95av2.db","igraph","limma","purrr",
#"RColorBrewer","RCy3","readxl","RedeR", "scales","stats",
#"stringr","tibble","tidyr","topGO","transcriptogramer,"RTN""

if (!require("BiocManager")){
  install.packages("BiocManager")
}

if (!requireNamespace("transcriptogramer")) {
  BiocManager::install("transcriptogramer")
}

if (!requireNamespace("RTN")) {
  BiocManager::install("RTN")
}

if (!require("biomaRt")){
  BiocManager::install("biomaRt")
}

if (!require("circlize")) {
  install.packages("circlize")
}

if (!require("classInt")) {
  install.packages("classInt")
}


if (!require("dplyr")) {
  install.packages("dplyr")
}

if (!require("edgeR")) {
  BiocManager::install("edgeR")
}

if (!require("factoextra")) {
  install.packages("factoextra")
}

if (!require("FactoMineR")) {
  install.packages("FactoMineR")
}

if (!require("ggplot2")) {
  install.packages("ggplot2")
}

if (!require("ggrepel")) {
  install.packages("ggrepel")
}

if (!require("grid")) {
  install.packages("grid")
}

if (!require("gridExtra")) {
  install.packages("gridExtra")
}

if (!require("igraph")) {
  install.packages("igraph")
}

if (!require("purrr")) {
  install.packages("purrr")
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
}

if (!require("RCy3")) {
  BiocManager::install("RCy3")
}

if (!require("RedeR")) {
  BiocManager::install("RedeR")
}

if (!require("refGenome")) {
  install.packages("refGenome")
}

if (!require("scales")) {
  install.packages("scales")
}

if (!require("stats")) {
  install.packages("stats")
}

if (!require("tibble")) {
  install.packages("tibble")
}

if (!require("tidyr")) {
  install.packages("tidyr")
}

if (!require("topGO")) {
  install.packages("topGO")
}

if (!require("wesanderson")) {
  install.packages("wesanderson")
}



#Create folder estructure
if (!dir.exists("figuras")){
  dir.create("./figuras")
}

if (!dir.exists("terms")){
  dir.create("./terms")
}

if (!dir.exists("samples")){
  dir.create("./samples")
}

if (!dir.exists("levels")){
  dir.create("./levels")
}

if (!dir.exists("log")){
  dir.create("./log")
}

if (!dir.exists("tmpGraf")){
  dir.create("./tmpGraf")
}
