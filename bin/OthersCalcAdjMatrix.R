rm(list = ls())

#Did you change it to your base location?
baseDir="~/LeadNPC/"
setwd(baseDir)
source(file = "bin/00base.R")

load("tabelasFig1b.RData")

matAd<-matrix(nrow = 11, ncol = 11)
i=2
j=3
for(i in 1:11){
  for (j in i:11) {
    soma<-sum(c_filter$Freq[(c_filter$Var1 == i &
                              c_filter$Var2 == j)|
                              (c_filter$Var1 == j &
                              c_filter$Var2 == i)])
    matAd[i,j]<-soma
    matAd[j,i]<-soma
  }
}
write.csv(x = matAd,file = "./figures/AdjMatrixT1.csv")
size<-t(unique(c_filter[,c(1,4)]))
write.csv(x = size,file = "./figures/AdjMatrixT1.csv")

rm(list = ls())

load("tabelasFig1d.RData")
c_filter<-tabela_t2

matAd<-matrix(nrow = 11, ncol = 11)
i=2
j=3
for(i in 1:11){
  for (j in i:11) {
    soma<-sum(c_filter$Freq[(c_filter$Var1 == i &
                               c_filter$Var2 == j)|
                              (c_filter$Var1 == j &
                                 c_filter$Var2 == i)])
    matAd[i,j]<-soma
    matAd[j,i]<-soma
  }
}
write.csv(x = matAd,file = "./figures/AdjMatrixT2.csv")
