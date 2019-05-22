setwd("/home/clovis/Dropbox/Chumbo/")

load(file = "./Data/associationHs700.RData")
temp1<-data.frame(c1=as.integer(substr(association[,1],10,25)),
                  c2=as.integer(substr(association[,2],10,25)),
                  stringsAsFactors = F)

temp2<-as.data.frame(t(apply(temp1,MARGIN = 1,FUN = function(x){
  #cat(x[1],x[2],"\n")
  if(x[1]>x[2]){
    return(c(x[2],x[1]))
  }else{
    return(c(x[1],x[2]))
  }
})))
nrow(unique(temp2))

assocNoDup<-association[duplicated(temp2),]

#checagem
temp1<-data.frame(c1=as.integer(substr(assocNoDup[,1],10,25)),
                  c2=as.integer(substr(assocNoDup[,2],10,25)),
                  stringsAsFactors = F)

temp2<-as.data.frame(t(apply(temp1,MARGIN = 1,FUN = function(x){
  #cat(x[1],x[2],"\n")
  if(x[1]>x[2]){
    return(c(x[2],x[1]))
  }else{
    return(c(x[1],x[2]))
  }
})))
nrow(unique(temp2))

#checagem 2
tmp<-unique(assocNoDup[,1])[1:100]
for(i in 1:100){
  cat(tmp[i],": col1:",
    length(grep(tmp[i], assocNoDup[,1])),
    " col2:",
    length(grep(tmp[i], assocNoDup[,2])),
    "\n")
}

save(assocNoDup,file = "./Data/assocNoDup.RData")

